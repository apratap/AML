# #source("http://bioconductor.org/biocLite.R")
# #biocLite("limma")
# 
# library(limma)
# 
# gpr_files_dir <- "/work/DAT_118__AML/Analysis/dset1//microRNA//data/"
# list.files(gpr_files_dir,full.names=T)
# ?readTargets
# ?read.maimages

#source("http://bioconductor.org/biocLite.R")
#biocLite("MiChip")
library("MiChip")
setwd("~/dev//AML/")
source("MiChip_mod_functions.R")
#rawData_testing <- parseRawData_mod(datadir="/work/DAT_118__AML//Analysis/dset1//microRNA//data/testing/", pat="gpr")

#parse the gpr files
#returns eSet object
rawData <- parseRawData_mod(datadir="/work/DAT_118__AML//Analysis/dset1//microRNA//data/", pat="gpr")

#Removing Unwanted Rows and Correcting for Flags
cleanedRawData <- removeUnwantedRows(rawData, c("blank", "Blank"))

#Removes all empty spots, control spots, U6 RNA, non human spots from an ExpressionSet in the standard fashion. 
#A wrapper for removeUnwantedRows
cleanedRawData <- standardRemoveRows(cleanedRawData)

#showMethods(class=class(cleanedRawData))
#flagcorrected_cleanedRawData <- correctForFlags(cleanedRawData, intensityCutoff = 50)
#flagcorrected_cleanedRawData

#Summarizing Intensities
summedData <- summarizeIntensitiesAsMedian(cleanedRawData, minSumlength=0, madAdjust=0)
miRNA_exp_matrix <- exprs(summedData)


#map back the col names to patientIDs
library(synapseClient)
synapseLogin()
query_df <- synQuery("select * from entity where parentId == 'syn2330288'")
colnames(query_df) <- sub('entity.', '',colnames(query_df))
#change the col names of the expression matrix r
#map the marcucci id to patient id
patientIds <- unlist(lapply(colnames(m), function(x){ patientID <- query_df[query_df$name == paste0(x,'.gpr'),]$PatientID}))
colnames(miRNA_exp_matrix) <- patientIds

#log2 intensity boxplot
boxplotData(miRNA_exp_matrix, "dset1_miRNA_expression", "raw_intesity")


#this gives the mapping to probeID(sequence) to NAME
#currently in the expression matrix I am using probeID/sequence as some probeID's have no name associated
head(pData(featureData(summedData)))


#create provenance and upload to synapse
syn_miRNA_files <- synapseQuery('SELECT id, name FROM entity WHERE parentId=="syn2330288"')

## LINK THESE TOGETHER WITH PROVENANCE
act  <- Activity(used=syn_mRNA_files$entity.id, executed=c('https://github.com/apratap/AML/blob/master/AML_dset1_mRNA_analysis.R'))
act <- synStore(act)
mRNA_arrayQC1 <- File("dset1_mRNA_array_intensities.png",synapseStore=T,parentId ="syn2329716")
mRNA_arrayQC1 <- synStore(mRNA_arrayQC1, activity=act)
