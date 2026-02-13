# Codon Usage Counter (eggNOG → per COG codon usage)

This script counts **codon usage per CDS (gene)** from per-genome FASTA files and links genes to **COG IDs** using eggNOG-mapper annotations. It can write results either:

1) **Per COG** (default): one TSV per COG containing all genes annotated with that COG from all genomes 
2) **Per genome**: one TSV per genome containing all COG-annotated genes in each genome

It supports two annotation input modes:

- **per_genome**: one `.emapper.annotations` file per genome  
- **combined**: a single combined annotation file covering all genomes

---

## Requirements

- Python 3.x
- No external Python packages required

---

## Input data layout

### 1) CDS FASTA files (required)
A folder containing nucleotide CDS FASTA files named:

```
GENOMEID.cds.fna
```

Each record header must start with a gene ID that matches the first column of eggNOG output:

```
>gene1234 some optional description
ATG...
```

Example folder:

```
cds_folder/
  Acido1020.cds.fna
  GenomeB.cds.fna
  ...
```

### 2) eggNOG-mapper annotations (required)

#### Mode A: per-genome annotations (`--mode per_genome`)
Provide a folder containing per-genome eggNOG files named:

```
GENOMEID.emapper.annotations
```

Example:

```
eggnog_folder/
  Acido1020.emapper.annotations
  GenomeB.emapper.annotations
  ...
```

#### Mode B: combined annotations (`--mode combined`)
Provide a single combined eggNOG annotations file:

```
combined.emapper.annotations
```

---

## What COG assignments are used?

The script reads eggNOG “COG” assignments from **column 5** of the `.emapper.annotations` file (tab-delimited).  
It keeps entries whose ID starts with `COG` and **strips any taxonomic level suffix** (anything after `@`).It will basically count only root level.

Example field:

```
COG0457@1|root,COG0457@2|Bacteria,4NE6G@976|Bacteroidetes
```

The script will keep only:

- `COG0457`

and ignore lineage-specific non-COG IDs like `4NE6G`.

If a gene maps to multiple COGs, it will be included once per COG (same codon counts, multiple rows).

---

## Usage

Make the script executable (optional):

```bash
chmod +x codon_counts.py
```

### A) Per-genome eggNOG files → output by COG (default)
```bash
./codon_counts.py cds_folder/ eggnog_folder/ outdir/ --mode per_genome
```

### B) Combined eggNOG file → output by COG (default)
```bash
./codon_counts.py cds_folder/ combined.emapper.annotations outdir/ --mode combined
```

### C) Output one file per genome
Add `--output_style by_genome`:

```bash
./codon_counts.py cds_folder/ eggnog_folder/ outdir/ --mode per_genome --output_style by_genome
```

---

## Output formats

### 1) Output style: `by_cog` (default)
Creates one TSV per COG:

```
outdir/
  COG0001.tsv
  COG0457.tsv
  ...
```

Each file contains:

- `genome_geneID` = `GENOMEID_geneID`
- 64 codon count columns (AAA, AAT, ...)

Header example:

```
genome_geneID	AAA	AAT	...	CCC
```

### 2) Output style: `by_genome`
Creates one TSV per genome:

```
outdir/
  Acido1020.tsv
  GenomeB.tsv
  ...
```

Each file contains:

- `genome_geneID` = `GENOMEID_geneID`
- `COG`
- 64 codon count columns

Header example:

```
genome_geneID	COG	AAA	AAT	...	CCC
```

---

## Notes / gotchas

- Codons are counted in-frame from position 0 in each CDS, stepping by 3.
- Any codon containing characters outside `A/T/G/C` is ignored.
- Genes with no COG annotation are skipped.
- If you use `--mode per_genome` and a genome has no matching `.emapper.annotations` file, it will be skipped with a warning.

---

## Contact / attribution

Script maintained by: (add your name/lab here)
