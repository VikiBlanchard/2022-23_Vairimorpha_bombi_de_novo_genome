# Create the conda environment for RepeatModeler and RepeatMasker
# The RepeatModeler package contains the RepeatMasker package
> conda create -c bioconda -n repeatmodeler repeatmodeler
# Activate the environment for RepeatModeler
> conda activate repeatmodeler
# Install the dependencies for RepeatModeler
> conda update -c conda-forge perl-file-which
# Download the NINJA package for large-scale neighbor-joining phylogeny inference and
clustering
> mkdir bin
> cd bin
> wget https://github.com/TravisWheelerLab/NINJA/archive/refs/tags/0.95-cluster_only.
tar.gz
> tar -zxvf 0.95-cluster_only.tar.gz
> cd NINJA-0.95-cluster_only/NINJA/
> make # Create the ‘‘Ninja’’ executable file
> pwd
# The pwd Linux command prints the current working directory path
# The standard output of the "pwd" command will be used as a parameter of RepeatModeler

# Repeat masking using known metazoan repeats with RepeatMasker:
#!/usr/bin/env bash
# You should use scaffold files as the input in RepeatMasker
for sample in ‘ls *fa‘;do
RepeatMasker -species metazoa -s -parallel 10 -xsmall -alignments $sample
done
# -s: sensitive
# -parallel 10: number of threads
# -xsmall: softmasking, that is, change the repeat regions into lowercase, rather than N

#  Identify previously unknown repeats in your genome assembly using RepeatModeler:
#!/usr/bin/env bash
#1. Create a Database for RepeatModeler
BuildDatabase -name CLR CLR_scaffold.fa
BuildDatabase -name ONT ONT_scaffold.fa
BuildDatabase -name Hifi Hifi_scaffold.fa
# -name: The name of the database to create
#2. Run RepeatModeler
RepeatModeler -database CLR -pa 10 -LTRStruct -ninja_dir /home/assembly/bin/NINJA-0.95-
cluster_only/NINJA
RepeatModeler -database ONT -pa 10 -LTRStruct -ninja_dir /home/assembly/bin/NINJA-0.95-
cluster_only/NINJA
RepeatModeler -database Hifi -pa 10 -LTRStruct -ninja_dir /home/assembly/bin/NINJA-0.95-
cluster_only/NINJA
# -database: prefix name of the database that is used in the BuildDatabase function
# -pa: number of threads
# -LTRStruct: runs the LTR structural discovery pipeline for discovering LTR retrotransposons
# -ninja_dir: specify the NINJA folder

#Repeat masking with RepeatMasker using species-specific repeats that were found by RepeatModeler
#!/usr/bin/env bash
RepeatMasker -lib CLR-families.fa -s -parallel 10 -xsmall -alignments CLR_scaffold.
fa.masked
RepeatMasker -lib ONT-families.fa -s -parallel 10 -xsmall -alignments ONT_scaffold.fa.
masked
RepeatMasker -lib Hifi-families.fa -s -parallel 10 -xsmall -alignments Hifi_scaffold.fa.
masked
# -lib: specify your species-specific repeat FASTA file produced by RepeatModeler
# -s: sensitive
# -xsmall: softmasking, that is, change the repeat regions into lowercase, rather than N

