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
from collections import defaultdict
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
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "gff3,genome"
        ]
        run(cmd)
        # Datasets writes "ncbi_dataset.zip" to CWD by default
        produced = Path("ncbi_dataset.zip")
        if not produced.exists():
            raise RuntimeError("Expected datasets to produce 'ncbi_dataset.zip' in the current directory.")
        produced.replace(zip_path)

    # Extract zip
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(extract_dir)

    return extract_dir

# ----------------------------
# Main pipeline
# ----------------------------

def main(file_a: Path, file_c: Path, out_base: Path):
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

    # 3) Process each genome
    missing = []
    for genome_short, genes in genome_to_genes.items():
        accession = normalize_accession(short_to_acc.get(genome_short, ""))
        if not accession:
            missing.append(genome_short)
            continue


        out_fasta = cds_dir / f"{genome_short}.cds.fna"
        if out_fasta.exists():
            # Skip if already done
            continue

        # Download dataset (gff + genome fasta)
        dataset_dir = download_ncbi_dataset(accession, downloads_dir)

        # Find genomic fasta (.fna) inside dataset
        fna_path = find_first_file(dataset_dir, suffixes=[".fna", ".fna.gz"])
        if fna_path is None:
            raise RuntimeError(f"Could not find a .fna file inside dataset for {accession} at {dataset_dir}")

        if fna_path.suffix == ".gz":
            raise RuntimeError(
                f"Found gzipped FASTA {fna_path} but this simple script does not gunzip.\n"
                "Easiest fix: run `gunzip -k <file>` or modify script to handle gz."
            )

        contigs = read_fasta_to_dict(fna_path)

        # Extract CDS sequences and write per-genome FASTA
        with out_fasta.open("w", encoding="utf-8") as out:
            for new_id, original in genes:
                try:
                    contig_id, start, end, strand = parse_prodigal_style(original)
                except Exception as e:
                    sys.stderr.write(f"[WARN] {genome_short}: could not parse '{new_id}': {e}\n")
                    continue

                seq = contigs.get(contig_id)
                if seq is None:
                    sys.stderr.write(
                        f"[WARN] {genome_short}: contig '{contig_id}' not found for gene {new_id}\n"
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

        print(f"[OK] {genome_short}: wrote {out_fasta}")

    if missing:
        sys.stderr.write(
            "\n[WARN] These genome short codes were present in File C but missing in File A:\n"
            + "\n".join(sorted(missing)) + "\n"
        )

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage:\n"
            "  python extract_cds_from_gtdb_accessions.py genome2abbrev.csv sequenceid_new2original.csv OUTPUT_DIR\n"
        )
        sys.exit(1)

    file_a = Path(sys.argv[1])
    file_c = Path(sys.argv[2])
    out_base = Path(sys.argv[3])

    main(file_a, file_c, out_base)
