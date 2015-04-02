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
source("~/dev/apRs/synapseHelpers.R")

library("synapseClient")
synapseLogin()


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


###Download data for OHSU plates
#1. OHSU plates run at FIMM
OHSU_plate_run_at_FIMM_synId <- "syn3458537"
OHSU_plate_run_at_FIMM <- temp_getData(OHSU_plate_run_at_FIMM_synId)

#2. OHSU plates run at OHSU
OHSU_plate_run_at_OHSU_synid <- "syn3458940"
OHSU_plate_run_at_OHSU <- temp_getData(OHSU_plate_run_at_OHSU_synid)

#3. OHSU plates run at Phoenix
OHSU_plate_run_at_Phoenix_synid <- "syn3458186"
OHSU_plate_run_at_Phoenix <- temp_getData(OHSU_plate_run_at_Phoenix_synid)


## Dose Response fit
OHSU_plates_rawData <- rbind(OHSU_plate_run_at_FIMM, OHSU_plate_run_at_OHSU, OHSU_plate_run_at_Phoenix)


#setting all the -ve normViab to 0
#OHSU_plates_rawData$normValue[OHSU_plates_rawData$normValue < 0] = 0

#fix cellLine names
OHSU_plates_rawData$cellLine[OHSU_plates_rawData$cellLine == 'Molm14'] = 'MOLM14'
OHSU_plates_rawData$cellLine = as.character(OHSU_plates_rawData$cellLine)

#1.
library("doMC")
registerDoMC(4)
OHSU_plates_doseResp_fit <- ddply(.data=OHSU_plates_rawData,
                                  .variable=c('plate_origin', 'run_at', 'plate_version','plateNum', 
                                              'run_by', 'run_date', 'set', 'cellLine' , 'assay',
                                              'plate_reader', 'read_time', 'Drug', 'replicate_plate', 
                                              'drug_replicate_across_plates'),
                                  .fun = function(x) {get_drugResponse_stats( x$Concentration, x$normValue, useLog=T) },
                                  .parallel = T
                                )

#upload doseRespone fit data for all OHSU plates
synUploadFile('OHSU_plates_doseResp_fit.tsv', 'syn3479857',
              used=c(OHSU_plate_run_at_FIMM_synId, OHSU_plate_run_at_OHSU_synid, OHSU_plate_run_at_Phoenix_synid),
              executed = '')


