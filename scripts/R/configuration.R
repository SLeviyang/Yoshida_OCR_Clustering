library(dplyr)
library(plyr)
library(Biostrings)
library(SummarizedExperiment)
library(glmnet)
library(JASPAR2018)
library(TFBSTools)
library(ggplot2)

HOMER_PATH <- "/Users/sivanleviyang/Dropbox/OSX_software/homer/bin"

source("data.R")
source("peak_clusters.R")
source("create_motif_scores_matrices.R")
source("utilities_HOMER.R")