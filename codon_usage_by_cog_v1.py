#!/usr/bin/env python3

import os
import sys
import argparse
from collections import defaultdict

# --------------------------------------------------
# 1. Generate all 64 codons
# --------------------------------------------------

def generate_codons():
    bases = ["A", "T", "G", "C"]
    return [a + b + c for a in bases for b in bases for c in bases]


CODONS = generate_codons()


# --------------------------------------------------
# 2. Count codons
# --------------------------------------------------

def count_codons(seq):
    seq = seq.upper().replace("\n", "").replace(" ", "")
    counts = dict.fromkeys(CODONS, 0)

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and all(base in "ATGC" for base in codon):
            counts[codon] += 1

    return counts


# --------------------------------------------------
# 3. Parse FASTA
# --------------------------------------------------

def parse_fasta(filepath):
    sequences = {}
    with open(filepath) as f:
        gene_id = None
        seq_lines = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if gene_id:
                    sequences[gene_id] = "".join(seq_lines)
                gene_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)

        if gene_id:
            sequences[gene_id] = "".join(seq_lines)

    return sequences


# --------------------------------------------------
# 4A. Parse per-genome eggNOG
# --------------------------------------------------

def parse_eggnog_file(filepath):
    gene_to_cogs = defaultdict(set)

    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            gene_id = parts[0]
            eggnog_field = parts[4]

            if eggnog_field == "-":
                continue

            entries = eggnog_field.split(",")

            for entry in entries:
                cog = entry.split("@")[0]
                if cog.startswith("COG"):
                    gene_to_cogs[gene_id].add(cog)

    return gene_to_cogs


# --------------------------------------------------
# 4B. Parse combined eggNOG file
# --------------------------------------------------

def parse_combined_eggnog(filepath):
    gene_to_cogs = defaultdict(set)

    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            gene_id = parts[0]
            eggnog_field = parts[4]

            if eggnog_field == "-":
                continue

            entries = eggnog_field.split(",")

            for entry in entries:
                cog = entry.split("@")[0]
                if cog.startswith("COG"):
                    gene_to_cogs[gene_id].add(cog)

    return gene_to_cogs


# --------------------------------------------------
# 5. Main
# --------------------------------------------------

def main(cds_folder, eggnog_input, output_folder, mode):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    cog_data = defaultdict(list)

    # If combined mode, parse once
    if mode == "combined":
        print("Parsing combined eggNOG file...")
        global_gene_to_cogs = parse_combined_eggnog(eggnog_input)

    # Iterate CDS files
    for filename in os.listdir(cds_folder):

        if not filename.endswith(".cds.fna"):
            continue

        genome_name = filename.rsplit(".cds.fna", 1)[0]
        cds_path = os.path.join(cds_folder, filename)

        print(f"Processing {genome_name}")

        sequences = parse_fasta(cds_path)

        # Choose eggNOG mapping source
        if mode == "per_genome":
            eggnog_file = os.path.join(
                eggnog_input,
                genome_name + ".emapper.annotations"
            )

            if not os.path.exists(eggnog_file):
                print(f"Warning: No eggNOG file for {genome_name}")
                continue

            gene_to_cogs = parse_eggnog_file(eggnog_file)

        else:  # combined
            gene_to_cogs = global_gene_to_cogs

        for gene_id, seq in sequences.items():

            if gene_id not in gene_to_cogs:
                continue

            counts = count_codons(seq)
            genome_gene_id = f"{genome_name}_{gene_id}"

            for cog in gene_to_cogs[gene_id]:
                row = [genome_gene_id] + [str(counts[c]) for c in CODONS]
                cog_data[cog].append(row)

    # Write output
    header = ["genome_geneID"] + CODONS

    for cog, rows in cog_data.items():
        outpath = os.path.join(output_folder, f"{cog}.tsv")
        with open(outpath, "w") as out:
            out.write("\t".join(header) + "\n")
            for row in rows:
                out.write("\t".join(row) + "\n")

    print("Done.")


# --------------------------------------------------
# 6. Argument parser
# --------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Count codon usage per gene and separate by COG."
    )

    parser.add_argument("cds_folder",
                        help="Folder with GENOMEID.cds.fna files")

    parser.add_argument("eggnog_input",
                        help="Either folder (per_genome mode) or combined eggNOG file")

    parser.add_argument("output_folder",
                        help="Output folder")

    parser.add_argument("--mode",
                        choices=["per_genome", "combined"],
                        required=True,
                        help="Use per_genome or combined eggNOG annotations")

    args = parser.parse_args()

    main(args.cds_folder,
         args.eggnog_input,
         args.output_folder,
         args.mode)
