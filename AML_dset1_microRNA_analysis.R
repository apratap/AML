#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")

library(limma)

gpr_files_dir <- "/work/DAT_118__AML/Analysis/dset1//microRNA//data/"

list.files(gpr_files_dir,full.names=T)

?readTargets

?read.maimages