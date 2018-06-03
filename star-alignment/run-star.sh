#! /bin/bash

# You run this file from the command line by passing in the name of the organism
# (fish, fly, or mouse), like so:
#
#    $ run-star.sh mouse

org="$1" # Stores the value of the organism name you passed in to the script

inpath="/data/reads/${org}"

for infn in $(ls ${inpath}/*.fastq.gz); do
  echo ""
  echo ""
  echo "======================================================================="
  echo $infn

  sid=$(basename "$infn" | cut -d '_' -f1)
  echo "  sampleid: $sid"

  outdir="/data/alignments/${org}/${sid}/"
  mkdir $outdir

  cmd="STAR \
    --genomeLoad LoadAndKeep \
    --genomeDir /data/star/${org} \
    --outFilterType BySJout \
    --quantMode GeneCounts \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --runThreadN 12 \
    --readFilesIn ${infn} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${outdir} \
    --outReadsUnmapped Fastx \
    --outSAMattributes NH HI AS nM MD XS \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMcompression 6 \
    --limitBAMsortRAM 15000000000"
  #echo $cmd
  eval $cmd
done

STAR --genomeLoad Remove --genomeDir /data/star/${org}

# Now index the files
for bam in $(ls /data/alignments/${org}/*/*.bam); do
  echo "Indexing: ${bam}"
  samtools index ${bam}
done
