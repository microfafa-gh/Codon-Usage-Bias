#!/usr/bin/env python3
"""
Extract CDS nucleotide sequences per genome using:
- File A: genome2abbrev.csv (maps short -> accession)
- File C: sequenceid_new2original.csv (maps short_gene -> prodigal-style coordinates)

Workflow:
1) For each genome short code found in File C:
   - Find its accession using File A
   - Download genome package (GFF + genomic FASTA) using NCBI Datasets CLI
   - Extract CDS sequences from the genomic FASTA based on coordinates in File C
   - Write per-genome CDS FASTA file
"""

import csv
import os
import re
import sys
import shutil
import zipfile
import subprocess
import argparse
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ----------------------------
# Small helper functions
# ----------------------------

def revcomp(seq: str) -> str:
    """Reverse-complement a DNA sequence (assumes ACGTN only; preserves N)."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def read_fasta_to_dict(fasta_path: Path) -> dict:
    """
    Read a FASTA file into a dict: {header_id: sequence}.
    header_id is the first token after '>' (up to first whitespace).
    """
    seqs = {}
    header = None
    chunks = []

    with fasta_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks)

    return seqs

def parse_prodigal_style(s: str):
    """
    Parse the "original" column example:
      NZ_JAEKFI010000001.1_2 # 499 # 2883 # 1 # ID=1_2;...

    Returns:
      contig_id (str), start (int), end (int), strand (int)

    Notes:
    - contig_id is obtained by removing the trailing "_<geneNumber>" from the first field
    - start/end are 1-based inclusive
    - strand is 1 or -1
    """
    parts = [p.strip() for p in s.split("#")]
    if len(parts) < 4:
        raise ValueError(f"Not enough '#' fields: {s}")

    first_field = parts[0]  # e.g. NZ_..._2
    start = int(parts[1])
    end = int(parts[2])
    strand = int(parts[3])

    # Remove trailing _<digits> to get contig name
    contig_id = re.sub(r"_(\d+)$", "", first_field)

    # Ensure start <= end
    if start > end:
        start, end = end, start

    if strand not in (1, -1):
        raise ValueError(f"Unexpected strand '{strand}' in: {s}")

    return contig_id, start, end, strand

def run(cmd, *, cwd=None):
    """Run a command, raise a nice error if it fails."""
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = (
            f"Command failed (exit {proc.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}\n"
        )
        raise RuntimeError(msg)
    return proc

def find_first_file(root: Path, suffixes):
    """Find the first file in 'root' matching any suffix in suffixes."""
    for suf in suffixes:
        hits = list(root.rglob(f"*{suf}"))
        if hits:
            return hits[0]
    return None

def normalize_accession(acc: str) -> str:
    """
    GTDB often prefixes accessions with RS_ or GB_.
    NCBI wants bare GCF_/GCA_ accessions.
    Examples:
      RS_GCF_016650635.1 -> GCF_016650635.1
      GB_GCA_123456789.1 -> GCA_123456789.1
    """
    acc = acc.strip()
    return re.sub(r"^(RS_|GB_)", "", acc)


def download_ncbi_dataset(accession: str, out_dir: Path) -> Path:
    """
    Download NCBI dataset zip for a genome accession using NCBI Datasets CLI.
    Includes GFF + genomic FASTA.

    Returns path to the extracted dataset directory.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = out_dir / f"{accession}.zip"
    extract_dir = out_dir / accession

    # Skip download/extract if already present
    if extract_dir.exists() and any(extract_dir.rglob("*.fna")):
        return extract_dir

    # Download zip
    if not zip_path.exists():
        tmp_cwd = Path(tempfile.mkdtemp(prefix=f"datasets_{accession}_", dir=str(out_dir)))
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "gff3,genome"
        ]
        try:
            run(cmd, cwd=tmp_cwd)
            # Datasets writes "ncbi_dataset.zip" to CWD by default
            produced = tmp_cwd / "ncbi_dataset.zip"
            if not produced.exists():
                raise RuntimeError("Expected datasets to produce 'ncbi_dataset.zip' in the temporary directory.")
            produced.replace(zip_path)
        finally:
            shutil.rmtree(tmp_cwd, ignore_errors=True)

    # Extract zip
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(extract_dir)

    return extract_dir


def detect_array_index() -> int:
    """Detect array task index from common HPC schedulers; return -1 if not found."""
    env_candidates = [
        "SLURM_ARRAY_TASK_ID",
        "PBS_ARRAY_INDEX",
        "SGE_TASK_ID",
        "LSB_JOBINDEX",
    ]
    for name in env_candidates:
        value = os.environ.get(name)
        if value is None:
            continue
        try:
            return int(value)
        except ValueError:
            continue
    return -1


def normalize_array_index(raw_index: int, array_total: int) -> int:
    """
    Accept either 0-based [0..array_total-1] or 1-based [1..array_total] index.
    Returns normalized 0-based index.
    """
    if 0 <= raw_index < array_total:
        return raw_index
    if 1 <= raw_index <= array_total:
        return raw_index - 1
    raise ValueError(
        f"array_index={raw_index} is invalid for array_total={array_total}. "
        f"Provide 0..{array_total - 1} or 1..{array_total}."
    )


def process_one_genome(genome_short: str, genes, accession: str, downloads_dir: Path, cds_dir: Path):
    """Download data for one genome and write its CDS FASTA."""
    out_fasta = cds_dir / f"{genome_short}.cds.fna"
    if out_fasta.exists():
        return "skipped", genome_short, 0, []

    warnings = []
    dataset_dir = download_ncbi_dataset(accession, downloads_dir)

    fna_path = find_first_file(dataset_dir, suffixes=[".fna", ".fna.gz"])
    if fna_path is None:
        raise RuntimeError(f"Could not find a .fna file inside dataset for {accession} at {dataset_dir}")

    if fna_path.suffix == ".gz":
        raise RuntimeError(
            f"Found gzipped FASTA {fna_path} but this script does not gunzip.\n"
            "Easiest fix: run `gunzip -k <file>` or modify script to handle gz."
        )

    contigs = read_fasta_to_dict(fna_path)
    written = 0

    with out_fasta.open("w", encoding="utf-8") as out:
        for new_id, original in genes:
            try:
                contig_id, start, end, strand = parse_prodigal_style(original)
            except Exception as e:
                warnings.append(f"[WARN] {genome_short}: could not parse '{new_id}': {e}")
                continue

            seq = contigs.get(contig_id)
            if seq is None:
                warnings.append(
                    f"[WARN] {genome_short}: contig '{contig_id}' not found for gene {new_id}"
                )
                continue

            # 1-based inclusive -> python slice
            subseq = seq[start - 1:end]
            if strand == -1:
                subseq = revcomp(subseq)

            header = f">{new_id} accession={accession} contig={contig_id} start={start} end={end} strand={strand}"
            out.write(header + "\n")

            # Wrap sequence to 60 chars/line
            for i in range(0, len(subseq), 60):
                out.write(subseq[i:i+60] + "\n")
            written += 1

    return "ok", genome_short, written, warnings

# ----------------------------
# Main pipeline
# ----------------------------

def main(
    file_a: Path,
    file_c: Path,
    out_base: Path,
    workers: int = 1,
    array_index: int = -1,
    array_total: int = 1,
):
    """
    file_a: genome2abbrev.csv
    file_c: sequenceid_new2original.csv
    out_base: output directory
    """
    out_base.mkdir(parents=True, exist_ok=True)
    downloads_dir = out_base / "downloads"
    cds_dir = out_base / "cds_fasta"
    cds_dir.mkdir(parents=True, exist_ok=True)

    # 1) Read File A: short -> accession
    short_to_acc = {}
    with file_a.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row["accession"].strip()
            short = row["short"].strip()
            if acc and short:
                short_to_acc[short] = acc

    if not short_to_acc:
        raise RuntimeError("No mappings found in File A. Check headers: accession,short,taxa")

    # 2) Read File C and group by genome short code
    # File C is shown as tab-separated in your head output; handle both TSV and CSV robustly.
    genome_to_genes = defaultdict(list)

    with file_c.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Try split by tab first; if not, fall back to comma; else whitespace split (last resort)
            if "\t" in line:
                cols = line.split("\t", 1)
            elif "," in line:
                cols = line.split(",", 1)
            else:
                cols = line.split(None, 1)

            if len(cols) < 2:
                continue

            new_id = cols[0].strip()      # e.g. Pseudo1_g2
            original = cols[1].strip()    # e.g. NZ_... # start # end # strand # ...

            # genome short code is prefix before first "_g"
            # (works for Pseudo1_g123 etc.)
            if "_g" not in new_id:
                continue
            genome_short = new_id.split("_g", 1)[0]

            genome_to_genes[genome_short].append((new_id, original))

    if not genome_to_genes:
        raise RuntimeError("No gene mappings found in File C. Check its delimiter and contents.")

    # 3) Join metadata and keep only genomes with accessions
    missing = []
    work_items = []
    for genome_short, genes in sorted(genome_to_genes.items()):
        accession = normalize_accession(short_to_acc.get(genome_short, ""))
        if not accession:
            missing.append(genome_short)
            continue
        work_items.append((genome_short, genes, accession))

    # Optional array sharding for cluster schedulers
    if array_total < 1:
        raise ValueError("--array-total must be >= 1")

    if array_total > 1:
        if array_index < 0:
            detected = detect_array_index()
            if detected < 0:
                raise ValueError(
                    "Array mode requested (--array-total > 1), but no array index provided. "
                    "Pass --array-index or set scheduler env var (SLURM_ARRAY_TASK_ID/PBS_ARRAY_INDEX/...)."
                )
            array_index = detected

        shard_idx = normalize_array_index(array_index, array_total)
        work_items = [
            item for i, item in enumerate(work_items)
            if i % array_total == shard_idx
        ]
        print(
            f"[INFO] Array shard {shard_idx}/{array_total}: processing {len(work_items)} genomes on this task."
        )

    if workers < 1:
        raise ValueError("--workers must be >= 1")

    # 4) Execute genome work (serial or parallel on one node)
    if workers == 1:
        for genome_short, genes, accession in work_items:
            status, gshort, written, warns = process_one_genome(
                genome_short, genes, accession, downloads_dir, cds_dir
            )
            for w in warns:
                sys.stderr.write(w + "\n")
            if status == "ok":
                print(f"[OK] {gshort}: wrote {cds_dir / f'{gshort}.cds.fna'} ({written} CDS)")
            else:
                print(f"[SKIP] {gshort}: {cds_dir / f'{gshort}.cds.fna'} already exists")
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            fut_to_short = {
                ex.submit(
                    process_one_genome,
                    genome_short,
                    genes,
                    accession,
                    downloads_dir,
                    cds_dir,
                ): genome_short
                for genome_short, genes, accession in work_items
            }
            for fut in as_completed(fut_to_short):
                status, gshort, written, warns = fut.result()
                for w in warns:
                    sys.stderr.write(w + "\n")
                if status == "ok":
                    print(f"[OK] {gshort}: wrote {cds_dir / f'{gshort}.cds.fna'} ({written} CDS)")
                else:
                    print(f"[SKIP] {gshort}: {cds_dir / f'{gshort}.cds.fna'} already exists")

    if missing:
        sys.stderr.write(
            "\n[WARN] These genome short codes were present in File C but missing in File A:\n"
            + "\n".join(sorted(missing)) + "\n"
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract CDS nucleotide FASTA files from GTDB accession mappings."
    )
    parser.add_argument("file_a", type=Path, help="genome2abbrev.csv")
    parser.add_argument("file_c", type=Path, help="sequenceid_new2original.csv")
    parser.add_argument("out_base", type=Path, help="Output directory")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers per node (default: 1).",
    )
    parser.add_argument(
        "--array-index",
        type=int,
        default=-1,
        help="Array task index (0-based or 1-based). Optional if scheduler env var exists.",
    )
    parser.add_argument(
        "--array-total",
        type=int,
        default=1,
        help="Total number of array shards/tasks (default: 1 = no sharding).",
    )
    args = parser.parse_args()

    main(
        args.file_a,
        args.file_c,
        args.out_base,
        workers=args.workers,
        array_index=args.array_index,
        array_total=args.array_total,
    )
