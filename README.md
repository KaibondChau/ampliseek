# AmpliSeek
This repository accompanies the manuscript Chau KK, Matlock W, Constantinides B, Lipworth S, Newbold L, Tipper H, Goodall T, Brett H, Hughes J, Crook DW, Eyre DW, Read DS, Walker AS, Stoesser N. Evaluation of metagenomic, 16S rRNA gene and ultra-plexed PCR-based sequencing approaches for profiling antimicrobial resistance gene and bacterial taxonomic composition of polymicrobial samples

AmpliSeek is a BBTools wrapper pipeline for mapping sequencing reads generated by the AmpliSeq for Illumina Antimicrobial Resistance Research Panel. AmpliSeek processes raw Illumina outputs files (fastq.gz) by trimming (BBMap), merging (BBMerge) and mapping stringently (BBMapSkimmer) to AmpliSeq target sequences. There is also an option to specify mapping thresholds manually. 


## Installation
Conda + pip
```
conda create -n ampliseek -c bioconda python=3 bbmap biopython
conda activate ampliseek
pip install git+https://github.com/KaibondChau/ampliseek
```

## Usage
| Required arguments | Description |
| --- | --- |
| `-f` or `--forward_reads` | Input forward fastq.gz file |
| `-r` or `--reverse_reads` | Input reverse fastq.gz file |

| Optional arguments | Description |
| --- | --- |
| `-h` | Displays help information |
| `-id` or `--id_filter` | Value to pass to BBMapSkimmer idfilter [0-1] |
