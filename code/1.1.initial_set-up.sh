#!/bin/bash 
## Usage: initial_set-up.sh isolate_name
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

#__appname__ = ''
#__author__ = 'Viki Blanchard (viki.blanchard.2018@live.rhul.ac.uk)'
#__version__ = '0.0.1'

#####################################
### Initialise output directories ###
#####################################

# Make directory for results
mkdir results_$1

# Make directory for reports 
mkdir reports_$1

# Make directory for log files 
mkdir reports_$1/logs

########################################
### Initialise anaconda environments ###
########################################

# 1.) Filter high quality reads
#conda create --yes -n read_quality_filter -c bioconda porechop NanoLyse NanoFilt NanoPlot 

# 3.) Canu assembly 
#conda create --yes -n canu_assembly -c bioconda canu quast

# 4.) BUSCO analysis 
#conda create --yes -n busco_analysis -c conda-forge -c bioconda busco=5.4.2 
