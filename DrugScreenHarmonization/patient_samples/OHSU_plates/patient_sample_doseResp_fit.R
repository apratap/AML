options(as.factor=F)
library("gdata")
library("reshape2")
library("tidyr")
library("plyr")
library("dplyr")
library("ggplot2")
library("GGally")

source("~/dev/apRs/drugScreen_functions.R")
source("~/dev/apRs/plotting_helpers.R")
source("~/dev/apRs/synapseHelpers.R")  #does the login
library("nplr")

tmp_iterator <- function(df){
  tryCatch({
    stats <- get_drugResponse_stats(df$Concentration, df$normValue, useLog=T)
  },error=function(e){
    print(dim(df))
    print(unique(df$cellLine))
    print(unique(df$Drug))
    print(unique(df$run_at))
    print(unique(df$assay))
    print(unique(df$plate_reader))
    print(unique(df$read_time))
    return(data.frame(nplr_error=TRUE))
    #stop(e)
  })
}

#remove unwanted row which have control , blank etc
temp_remove_rows <- function(df){
  filter(df, ! Drug %in% c('', 'nan', 'control', 'Control'))
}

remove_CTxG <- function(df){
  filter(df, assay != 'CTxG')
}

temp_getData <- function(synid){
  x <- synGet(synid)
  x <- read.csv(x@filePath, sep="\t", header=T)
  x <- temp_remove_rows(x)
  x <- remove_CTxG(x)
}

###Download raw compiled data for patient sample run on OHSU plates
run_at_phoenix_synid <- "syn4932400"
run_at_phoenix <- temp_getData(run_at_phoenix_synid)

# QC : testing the #replicates for QC purposes
drug_replicates <- run_at_phoenix %>%
  group_by(Drug, plate_version, run_at,run_date,
           assay, plate_reader, read_time, plateNum, set,drug_replicate_across_plates) %>%
  summarise(replicate_count = length(normValue))


#1.
library("doMC")
registerDoMC(12)
library(nplr)
doseResp_fit <- ddply(.data=run_at_phoenix,
                      .variable=c('Drug', 'plate_version', 'run_at', 'run_date',
                                  'assay', 'plate_reader', 'read_time', 'plateNum', 'set'),
                                      #all the replicates on a single plate are clubbed together
                                      .fun = function(x) {tmp_iterator(x) },
                                      .parallel = T )

outfile <- "OHSU_patient_samples_on_OHSU_plates_doseResp_fit.tsv"
write.table(doseResp_fit, file=outfile,sep="\t")
#upload doseRespone fit data for all OHSU plates
synStore(File(outfile, parentId = 'syn4932396'),
         used=c(run_at_phoenix_synid),
         executed = '')
unlink(outfile)
