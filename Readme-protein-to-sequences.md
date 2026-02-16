reads File A (genome2abbrev.csv) to map short -> accession
reads File C (sequenceid_new2original.csv) to map each short genome to its gene coordinates
downloads each genome package (using the NCBI Datasets CLI, simplest reliable way from an accession) including GFF + genomic FASTA (.fna)
extracts CDS nucleotide sequences from the genome FASTA using the Prodigal-style coordinates in File C
writes one FASTA per genome containing the CDS sequences
Requirements (one-time)
Install NCBI Datasets CLI (recommended)
This script expects the command datasets to be available.
On many systems you can install via conda:
```
conda install -c conda-forge ncbi-datasets-cli
```
(If you already have datasets, you’re good.)
To run :

python extract_cds_from_gtdb_accessions.py \
  genome2abbrev.csv \
  sequenceid_new2original.csv \
  out_cds

Parallel on one node:

python extract_cds_from_gtdb_accessions.py \
  genome2abbrev.csv \
  sequenceid_new2original.csv \
  out_cds \
  --workers 8

Cluster array mode (recommended for many genomes):
- Use `--array-total N` to split genomes into `N` shards.
- Each array task processes only its own shard.
- The script auto-detects task index from `SLURM_ARRAY_TASK_ID`, `PBS_ARRAY_INDEX`, `SGE_TASK_ID`, or `LSB_JOBINDEX`.

Example (SLURM):
```
#!/bin/bash
#SBATCH --job-name=cds_extract
#SBATCH --cpus-per-task=4
#SBATCH --array=0-19

python extract_cds_from_gtdb_accessions.py \
  genome2abbrev.csv \
  sequenceid_new2original.csv \
  out_cds \
  --workers 4 \
  --array-total 20
```

Example (PBS):
```
#PBS -J 1-20

python extract_cds_from_gtdb_accessions.py \
  genome2abbrev.csv \
  sequenceid_new2original.csv \
  out_cds \
  --workers 4 \
  --array-total 20
```

Notes:
- `--workers` = processes on a single node/task.
- `--array-total` = number of cluster shards/tasks.
- You can combine both (`array` across nodes + `workers` per node).
- If a per-genome output file already exists, that genome is skipped.
