#!/bin/bash 
## Usage: gene_cluster_summary_parser 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

##########################
### Read in data files ###
##########################

# Read in gene cluster summaries
gene_cluster_summary=${1}
gene_cluster_summary='/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/10.SynIma_runs/All_species_vs_1x/GENE_CLUSTERS_SUMMARIES.OMCL/GENE_CLUSTERS_SUMMARIES.clusters_and_uniques'

# Read in gff3 
V_bombi_gff3=${2}
V_bombi_gff3='/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/7.SynIma_files_Vairimorpha_bombi_8.1-3/Vairimorpha_bombi_8.1-3no1x/Vairimorpha_bombi_8.1-3no1x.annotation.gff3'

# Read in genome.fasta
V_bombi_genome=${3}
V_bombi_genome='/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/7.SynIma_files_Vairimorpha_bombi_8.1-3/Vairimorpha_bombi_8.1-3no1x/Vairimorpha_bombi_8.1-3no1x.genome.fa'


################################################
### 1.) Save all contig IDs & gene IDs to df ### 
################################################

# 

############################################################
### 3.) Save V.bombi genes with orthologs and contig IDs ###
############################################################

# Remove singletons and empty lines
Tidy_cluster_frame=`awk '$0 !~ /uniq_/' ${gene_cluster_summary} | awk NF`

# Save orthogroups with V.bombi and other species 
Tidy_IDs=`awk '$0 !~ /uniq_/' ${gene_cluster_summary} | awk NF | awk '{print $1}' | sort | uniq `

# Remove self matches 
for each_orthogroup in ${Tidy_IDs}; do 

    # Save number of species in orthogroup
    number_of_species_in_orthogroup=`awk -v ortho_ID=${each_orthogroup} '$1 == ortho_ID' ${gene_cluster_summary} | awk '{print $2}' | sort | uniq | wc -l`

    # Save orthogroup IDs without self-matches
    if [[ $number_of_species_in_orthogroup > 1 ]]; 
    then
        echo "$each_orthogroup"; 
    fi

done 




# Get orthogroups and gene IDs for V.bombi only
awk '$2 ~ /Vairimorpha_bombi_8.1-3no1x/ {print $1 "\t" $4}' > Vairi_ortho_summary.tab




# Save orthogroups with V.bombi and other species 
orthogroup_ids=`awk '{print $1}' Vairi_ortho_summary.tab | sort | uniq `

# Save number of species per orthogroup
for each_orthogroup in ${orthogroup_ids}; do 
    number_of_species_in_orthogroup=`awk -v ortho_ID=${each_orthogroup} '$1 == ortho_ID' ${gene_cluster_summary} | sort | uniq | wc -l`
    awk -v ${number_of_species_in_orthogroup}
done 

awk -v username='some username' -v line=32 'NR == line { $0 = $0 username } 1' file

awk '$0 !~ /uniq_/' ${gene_cluster_summary} | awk '$2 ~ /Vairimorpha_bombi_8.1-3no1x/ ' 

# Count unique species + header
awk '(index($0) {print})' ${gene_cluster_summary}
awk '{print $2}' ${gene_cluster_summary} | sort | uniq | wc -l

# Skip lines with unique orthogroups

    # Check if orthogroup contains V.bombi 
    #elif [[ $line == *"Vairimorpha_bombi_8.1-3no1x"]] then 
    #    awk 'NR > 1{print$1}'
    #fi    

# Save contig IDs for orthologous genes 

# Save summary table

# Look at GC content of every contig
seqkit fx2tab --name --only-id --gc 'results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta' | cat > "results_Vairimorpha_bombi_8.1-3/${isolate}_content_per_tig.tsv"