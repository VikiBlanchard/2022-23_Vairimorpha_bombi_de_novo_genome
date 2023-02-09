#!/bin/bash 
## Usage: 4.1.0.run_busco_analysis.sh 
## Inputs: isolate_name, PATH/TO/new_assembly.fasta
## Outputs: 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

#######################################################
### Set up conda environment to run the analysis in ###
#######################################################

# Activate the environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate busco_analysis

# Save path to new assembly into variable 
#path_to_assembly=${path_to_no_bombus_assembly}
path_to_assembly=${2}
isolate=${1}
#path_to_assembly='/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8_1/3.no_bombus_Vairimorpha_bombi_8_1_assembly/test-Vb.contigs.fasta'

########################################
### Run BUSCO analysis on new genome ###
########################################

# Make output directory 
mkdir results_${isolate}/4.busco_analysis_${3}
step_4_results_dir="results_${isolate}/4.busco_analysis_${3}"

# Move to logs directory in reports
mkdir reports_${isolate}/buscologs_${3}
cd reports_${isolate}/buscologs_${3}

# Run BUSCO analysis on new assembly for microsporidia database
busco -m genome -i ${path_to_assembly} -l microsporidia_odb10 -o ../../${step_4_results_dir}/"microsporidia"

# Run BUSCO analysis on new assembly for fungal database
busco -m genome -i ${path_to_assembly} -l fungi_odb10 -o ../../${step_4_results_dir}/"fungi"

##################################
### Make BUSCO summary reports ###
##################################

# Plot busco summary charts and send to reports directory
generate_plot.py -wd ../../$step_4_results_dir/"microsporidia"
mv ../../$step_4_results_dir/"microsporidia"/"busco_figure.png" ../"microsporidia_busco_figure.png"

generate_plot.py -wd ../../$step_4_results_dir/fungi/run_fungi_odb10

generate_plot.py -wd '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8/4.busco_analysis_results/fungi/run_fungi_odb10'
mv ../../$step_4_results_dir/"fungi"/"busco_figure.png" ../"fungi_busco_figure.png"
# Send summary plot to reports directory


# Send summary plot R file to code directory
mv ../../$step_4_results_dir/busco_analysis_results/busco_figure.R ../../code/4.1.1.busco_figure_${isolate}_${3}.R

# Correct R plot output path
sed -i "s/\"..\/..\/results_*\/4\.isolated_microsporidia_sequences\/busco_analysis\_results/my_output <- paste(\"..\/reports\_${isolate}/g" ../../code/4.1.1.busco_figure_${isolate}.R

# Move back to base working directory
cd ../..

# Deactivate the conda environment
conda deactivate 

