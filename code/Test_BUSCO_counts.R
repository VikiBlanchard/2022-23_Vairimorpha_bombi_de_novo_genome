#!/usr/bin/env Rscript
#' @title Reading and Parsing Data in R
#' @description Example code for reading in and manipulating data
#' @param 
#' @author Viki Webster

##########################
### Set up environment ###
##########################

### Clear global environment 
rm(list = ls())

#install.packages("rstudioapi")
library ("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
