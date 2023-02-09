Order to run Scripts on raw Nanopore data to make a de novo assembly

##Step 1##
In all scripts change ISOLATE=Name of isolate

##Step 2##
Make a .txt file called "${ISOLATE}_fastq_list.txt" containing all of the names of the fastq files that are in the your rawdata depositary:

/rds/general/user/shemming/projects/fisher-aspergillus-rawdata/live/minION/${ISOLATE}

##Step 3##
Run
qsub highqual_reads_filter.sh

Check nanoplot graphs made at the end

##Step 4##
run after changing any perameters you name (e.g. isolate name, genome size) 
qsub canu_assembly.sh 

##Step 5##
run 
qsub run_quast.sh

To find view graphs showing quality of assembly
