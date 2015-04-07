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



#1. get the FIMMS FO3D data for K562 cellLine
FIMMS_FO3D_K562 <- "syn3354588"
k562_FO3D <- synGet(FIMMS_FO3D_K562)
k562_FO3D <- read.xls(k562_FO3D@filePath)
#delete columns not needed
k562_FO3D$X <- NULL


#mod column names
colnames(k562_FO3D) <- gsub('_rawIntensity', '', colnames(k562_FO3D))
colnames(k562_FO3D) <- gsub('_..._assay.*_D.','',colnames(k562_FO3D), perl=T)
colnames(k562_FO3D) <- gsub('H1FO3D.D.','',colnames(k562_FO3D), perl=T)
colnames(k562_FO3D) <- gsub('\\.','_',colnames(k562_FO3D), perl=T)

k562_FO3D <- melt(k562_FO3D, id.vars=c('DWell', 'ProductId', 'ProductName', 'conc__num_',
                                'content', 'row', 'plate_layout'))
k562_FO3D <- k562_FO3D %>%
                separate(variable , into = c('center', 'cellLine', 'drugPanel', 
                                             'set_number', 'assay',
                                             'dose', 'plateNum'), 
                         sep="_", remove=T)


### MOLM CellLine
# get the FIMMS FO3D data for MOLM cellline
FIMMS_FO3D_MOLM <- "syn3375418"
molm_FO3D <- synGet(FIMMS_FO3D_MOLM)
molm_FO3D <- read.xls(molm_FO3D@filePath)

#mod column names
colnames(molm_FO3D) <- gsub('_rawIntensity', '', colnames(molm_FO3D))
colnames(molm_FO3D) <- gsub('H1FO3D.D.','',colnames(molm_FO3D), perl=T)
colnames(molm_FO3D) <- gsub('\\.','_',colnames(molm_FO3D), perl=T)
colnames(molm_FO3D) <- gsub('assay.*D.', '', colnames(molm_FO3D), perl=T)
colnames(molm_FO3D) <- gsub('_...__', '_',colnames(molm_FO3D), perl=T)
molm_FO3D$X <- NULL


molm_FO3D <- melt(molm_FO3D, id.vars=c('DWell', 'ProductId', 'ProductName', 'conc__num_',
                                       'content', 'row', 'plate_layout'))
molm_FO3D <- molm_FO3D %>%
  separate(variable, into = c('center', 'cellLine', 'num', 'drugPanel', 
                              'set_number', 'assay','dose', 'plateNum'), 
           sep="_", remove=T)

#cellLine + num is one column
molm_FO3D$cellLine <- paste0(molm_FO3D$cellLine,molm_FO3D$num)
molm_FO3D$num <- NULL


# combined raw data
FIMM_plates_rawData <- rbind(k562_FO3D, molm_FO3D)
FIMM_plates_rawData['plate_origin'] = 'FIMM'
FIMM_plates_rawData['run_at'] = FIMM_plates_rawData$center
FIMM_plates_rawData$center <- NULL

#1. calc norm Factors /plate
# taking the median of all the pos and neg wells
temp_calc_normFactors <- function(df){
  if(nrow(df) != 384){
    stop('the df passed doesnt have 384 rows .. check if this represents a full plate')
  }
  pos_control = median(df[df$content=="pos",'value'])
  neg_control = median(df[df$content=="neg",'value'])
  res = data.frame('pos_control'=pos_control,
                   'neg_control'=neg_control)
  return(res)
}

plateNormFactors  <- ddply(.data=FIMM_plates_rawData,
                           .variables=c('plateNum', 'assay', 'dose', 'cellLine'),
                           .fun = temp_calc_normFactors)


#merge norm 
FIMM_plates_rawData <- merge(FIMM_plates_rawData, plateNormFactors)
#create norm value
FIMM_plates_rawData['normValue'] = (FIMM_plates_rawData$value - FIMM_plates_rawData$pos_control) / (FIMM_plates_rawData$neg_control - FIMM_plates_rawData$pos_control)
  
#create conc column
FIMM_plates_rawData['conc'] <- apply(FIMM_plates_rawData, 1, function(x){
                                    factor <- 10^( 5 - as.numeric(gsub('D','',x['dose'])))
                                    as.numeric(x['conc__num_'])/factor
                                    })

#convert the nano molar conc to micro molar conc. (this matches to OHSU plates conc)
FIMM_plates_rawData$conc <- FIMM_plates_rawData$conc / 1000



outfile <- 'FIMMS_plates_rawData.tsv'
write.table(FIMM_plates_rawData, file=outfile, sep="\t", col.names=T, 
            quote = F, row.names=F)
FIMMS_cleaned_rawData_file <- File(outfile, parentId = 'syn3354585') 
FIMMS_cleaned_rawData_file <- synStore(FIMMS_cleaned_rawData_file,
                                       used = c('https://ohsu.box.com/s/93yuhq286g2i3280fa7wa1a4vf2us79w',
                                                'https://ohsu.box.com/s/p517k71bnq238w1n5ufx1hyo4st0ndw6'),
                                       executed = )
unlink(outfile)






#2. dose resp curve fit
library(doMC)
registerDoMC(4)
FIMM_plates_doseResponse_fit <- ddply(.data=filter(FIMM_plates_rawData, ProductId != ''),
                                      .variable=c('ProductId', 'ProductName', 'cellLine', 'plate_origin',
                                                  'run_at', 'assay','set_number'),
                                      .fun = function(x) { get_drugResponse_stats(x$conc, x$normValue, useLog=T) },
                                      .parallel = T
                                     )

#2. dose - response fit data
outfile <- 'FIMM_plates_doseResponse_fit.tsv'
write.table(doseResponse_fit_data, file=outfile, sep="\t", col.names=T, quote=F, row.names=F)
FIMM_plates_doseResponse_fit_file <- File(outfile, parentId = 'syn3375637')
FIMM_plates_doseResponse_fit_file <- synStore(FIMM_plates_doseResponse_fit_file,
                                              used = c(FIMMS_FO3D_K562, FIMMS_FO3D_MOLM))
unlink(outfile)




