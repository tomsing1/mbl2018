# Overview

A number of FASTQ files were generated beforehand to grease the skids
for data generation and analysis during the workshop.

This document outlines how these data were aligned to the genome and quantitated
at the gene-level using the [STAR aligner](https://github.com/alexdobin/STAR).
For ease/expediency, you will only be using
[Kallisto](https://pachterlab.github.io/kallisto/about) for transcript
quantitation during the workshop, but we wanted to show you the differences
between these two approaches for your future reference.

Note that there was a good deal of massaging the data before hand to make the
downstream processing tasks easier. These data have been places into the
`s3://mbl.data` Amazon S3 bucket.

# Data Overview

For this exercise, the directories of interest within the `s3://mbl.data`
Amazon S3 bucket are
the following:

```
s3://mbl.data/
├── references
│   ├── fish/
│   │   ├── gene_annotation.gtf.gz
│   │   └── genome.fa.gz
│   ├── fly/
│   │   ├── gene_annotation.gtf.gz
│   │   └── genome.fa.gz
│   ├── mouse/
│   │   ├── gene_annotation.gtf.gz
│   ┴   └── genome.fa.gz
├── reads
│   ├── may/
│   │   ├── fish/
│   │   │   ├── drdgwt01_R1.fastq.gz 
│   │   │   ├── drdgwt02_R1.fastq.gz 
│   │   │   ├── drdgwt03_R1.fastq.gz 
│   │   │   └── ...
│   │   ├── fly/
│   │   │   ├── dmneaohi01_R1.fastq.gz
│   │   │   ├── dmneaohi02_R1.fastq.gz 
│   │   │   ├── dmneaohi03_R1.fastq.gz
│   │   │   └── ...
│   │   ├── mouse/
│   │   │   ├── mmchko02_R1.fastq.gz
│   │   │   ├── mmchko03_R1.fastq.gz 
│   │   │   ├── mmchko04_R1.fastq.gz 
┴   ┴   ┴   └── ...
```

# Amazon Machine Setup

We first spun up a rather beefy Amazon machine, using the AMI prepared for
use in this workshop.

Spinning up a large instance to run files through STAR:

    AMI: mbl_2018_v1_ref - ami-260d7a59
    type: c5d.4xlarge (vCPU 16; RAM 32gb)

Worfklow:

1. Setup machine after launch

```bash
ssh -i ~/.ssh/mblsteve.pem.txt ubuntu@34.205.8.110
curl https://raw.githubusercontent.com/lianos/dotfiles/master/screenrc > .screenrc
screen -RD

## Install STAR
conda install star

## Installing samtools normally doesn't work, need to do this:
# https://github.com/bioconda/bioconda-recipes/issues/8210
conda create -c conda-forge -c bioconda -c default -n samtooltest samtools ncurses

## Setup data and result directories -- use the SSD
sudo mkdir /data
sudo mkfs -t ext4 /dev/nvme1n1
sudo mount /dev/nvme1n1 /data
sudo chown -Rf ubuntu /data

# Make directories to hold the data and the alignments -------------------------
# hold genome fastq files in /data/genomes/* directories
mkdir -p /data/genomes/mouse /data/genomes/fly /data/genomes/fish

# star genome indices of the genome fastqs will go in /data/star/*
mkdir -p /data/star/mouse /data/star/fly /data/star/fish

# We will download th fastq files by organism into /data/reads/*
mkdir -p /data/reads/mouse /data/reads/fly /data/reads/fish

# star read alignments will go in animal subdirectories
mkdir -p /data/alignments/mouse /data/alignments/fly /data/alignments/fish
```

# Process files per organism

```bash
cd /data
ORG="mouse"

# Create STAR index for organism ===============================================
# 1. Copy required annotatoion fiesgenome files --------------------------------
aws s3 cp \
  s3://mbl.data/references/${ORG}/genome.fa.gz \
  genomes/${ORG}/
gunzip genomes/${ORG}/genome.fa.gz

aws s3 cp \
  s3://mbl.data/references/${ORG}/gene_annotation.gtf.gz \
  genomes/${ORG}/
gunzip genomes/${ORG}/gene_annotation.gtf.gz

# Create STAR index ------------------------------------------------------------
mkdir star/${ORG}

STAR --runThreadN 8 --runMode genomeGenerate \
  --sjdbOverhang 50 --genomeDir star/${ORG} \
  --genomeFastaFiles genomes/${ORG}/genome.fa \
  --sjdbGTFfile genomes/${ORG}/gene_annotation.gtf
```

## Process Experimental Data ===================================================

We now have to download all of the data for this organism to the local machine.
We can use this using the `aws s3 sync` command to grab all of the files from
a locaiton in an Amazon S3 bucket and copy them directly to our machine.

The following command will download the directory that holds the
${ORG}-anism specific data files into a our local `/data/reads/${ORG}`
directory.

```bash
aws s3 sync s3://mbl.data/reads/may/${ORG} /data/reads/${ORG}
```

When the above command completes, we will have all of the FASTQ files for
this particular organism in our local /data/reads/${ORG} folder.

We wrote a `run-star.sh {$ORG}` bash script that will iterate through each
of the FASTQ files in `/data/reads/${ORG}` and align them with STAR:

```bash
/data/run-star.sh $ORG
```

You will now have all of the STAR output results per organism here in the
`/data/alignments/${ORG}/*` directories.

We can now remove the raw FASTQ files we downloaded.

```bash
rm -Rf /data/reads/$ORG
```

# Finally

When we run the above stes for each organism, ie.`$ORG="fly"`, `$ORG="fish"`,
and `$ORG="mouse"`, we can sync the results back up to the S3 bucket:

```bash
aws s3 sync alignments s3://mbl.data/star-alignments/may
```

