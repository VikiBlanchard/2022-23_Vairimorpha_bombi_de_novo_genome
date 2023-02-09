#!/bin/bash -l

### rfarrer@broadinstitute.org

# Exit with message
die() { echo "$@" ; exit 1; }

# Check input files are specified
if [ $# != 3 ] ; then
	echo "Usage: $0 <reference sequence> <Reads in FASTQ> <Reads2 or 'unpaired'>"
	die "       notes: If reference sequence hasn't been indexed, any files called the same name with a '.any suffix' will be deleted prior to indexing."
fi

# Check input files are readable
[ -r $1 ] || die "File: $1 does not appear to be valid"
[ -r $2 ] || die "File: $2 does not appear to be valid"
[ $3 == "unpaired" ] || [ -r $3 ] || die "Reads3 not set as 'unpaired' or readable"
SAMSCRIPT="/code/2.1.1.SAM_how_many_reads_align.pl"

# Output files
SAI="$2.sai"
SAM="$2-mem.sam"
BAM="$2-mem.bam"
SORTED="$2-mem.sorted"
SORTEDBAM="${SORTED}.bam"
SORTEDBAM="${2}-mem.sorted.bam"
CLEAN="${R1}-with-pairs-mem.sorted_clean.bam"

# Index the reference FASTA. First, check if it has already been indexed.
if [ ! -r $1.ann ] && [ ! -r $1.bwt ]  ; then
	echo "Indexing reference and first removing any previous associated file..."
	rm $1.*
	bwa index -a is $1
fi

# Align reads
echo "Making alignments with short reads using BWA mem..."
if [ $3 == "unpaired" ] ; then
	bwa mem $1 $2 -M > $SAM
else
	bwa mem $1 $2 $3 -M > $SAM
fi

# Count reads/nt aligned
#perl $SAMSCRIPT -s $SAM
				
# Convert SAM to BAM. Sort and make Pileup
echo "View, Sort, Index, Mpileup, and scripts..."
samtools view -S -b -u $SAM > $BAM
samtools sort $BAM -o $SORTEDBAM
samtools index $SORTEDBAM
samtools view -b -h -f 0x2 $SORTEDBAM > $CLEAN

# Cleanup
#rm $SAM
#rm $BAM	
