library("beadarray")
library("BeadArrayUseCases")
library("illuminaHumanv4.db")
library("synapseClient")

setwd("~/dev//AML")
dataDir <- "/work/DAT_118__AML/Analysis/dset1//mRNA/"


###get the Signal to Noise Ratio for ALL mRNA arrays
metrics <- read.table("/external-data/DAT_118__AML//dataset1-Raddich-lab//CARDINAL//mRNA/Metrics.txt",
                       sep="\t", header=F, as.is=T)
columns <- c("Date", "Beadchip",	"Section",	"FocusGrn",	"RegGrn",	"P05Grn",	"P95Grn",
             "FocusRed",	"RegRed",	"P05Red",	"P95Red")
colnames(metrics) <- columns
#get beadchip number
mRNA_beadlevelFiles <- list.files(dataDir, pattern="*.txt")
uniq_beadChips <- as.numeric(unique(gsub('^.*_(.+?)_.*$',"\\1",mRNA_beadlevelFiles)))

# not working as all the arrays dont have a value : self calculating
mRNA_metrics <- metrics[metrics$Beadchip %in% uniq_beadChips, ]
plot(mRNA_metrics$P95Grn / mRNA_metrics$P05Grn)




#read the beadlevel data
beadlevel_data <- readIllumina(dir=dataDir, useImages=F, illuminaAnnotation="Humanv4")
#illumina annotation based on the output of following
suggestAnnotation(beadlevel_data, verbose=T)


# testing
#beadlevel_data@sectionData
#class(beadlevel_data)
# slotNames(beadlevel_data)
#sample data
#beadlevel_data[[1]][1:10,]


#get the signal to noise ratio
mRNA_array_SNR <- lapply(beadlevel_data@beadData, function(x) { intensity <- (x[[4]]$Grn); r <- quantile(intensity,probs=c(.05,.95)); r['95%']/r['5%']})
png("dset1_mRNA_array_singal_to_noise_ratio.png", width=9,height=7, units="in", res=150)
plot(seq(1:length(mRNA_array_SNR)), mRNA_array_SNR, main="Signal to Noise ratio (AML dset1 mRNA Illumina beadarray)", 
     ylab='Ratio (95% / 5%) intensity', xlab="array index")
dev.off()


#boxplot of intensities across all the arrays
png("dset1_mRNA_array_intensities.png", width=9,height=7, units="in", res=150)
boxplot(beadlevel_data, transFun = logGreenChannelTransform, col = "green",
        ylab = expression(log[2](intensity)), las=2, outline=FALSE, main= "Array intensities")
dev.off()

#exressionQCPipeline
#produces QC for each array, time consuming and something that we might not need
#expressionQCPipeline(beadlevel_data,qcDir="/home/apratap/projects/AML/work/Analysis/dset1/mRNA/QC")


outlierplot(beadlevel_data, array=1)

#using the control information
combinedControlPlot(beadlevel_data)

#summarize
BSData <- summarize(beadlevel_data)
det = calculateDetection(BSData)

#add annotation
BSData <- addFeatureData(BSData, toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE","GENOMICLOCATION"))


png("dset1_mRNA_probe_detection_pvals.png", width=9,height=7, units="in", res=150)
boxplot(det, main='dset1_mRNA_probe_detection_pvals (Null: there is no expression) ', ylab='p-value(detection score)',xaxt='n', xlab='individual samples/array' )
dev.off()

###########
# Push THE PLOTS to synapse
###########
synapseLogin()
syn_mRNA_files <- synapseQuery('SELECT id, name FROM entity WHERE parentId=="syn2354333"')

## LINK THESE TOGETHER WITH PROVENANCE
act  <- Activity(used=syn_mRNA_files$entity.id, executed=c('https://github.com/apratap/AML/blob/master/AML_dset1_mRNA_analysis.R'))
act <- synStore(act)
mRNA_arrayQC1 <- File("dset1_mRNA_array_intensities.png",synapseStore=T,parentId ="syn2354332")
mRNA_arrayQC1 <- synStore(mRNA_arrayQC1, activity=act)

mRNA_arrayQC2 <- File("dset1_mRNA_probe_detection_pvals.png",synapseStore=T,parentId ="syn2354332")
mRNA_arrayQC2 <- synStore(mRNA_arrayQC2, activity=act)

mRNA_arrayQC3 <- File("dset1_mRNA_array_singal_to_noise_ratio.png",synapseStore=T,parentId ="syn2354332")
mRNA_arrayQC3 <- synStore(mRNA_arrayQC3, activity=act)


#log transformed expression vals
expVals <- exprs(BSData)
newColNames <- lapply(colnames(expVals), function(x){strsplit(x, '_')[[1]][1]})
colnames(expVals) <- newColNames       
raw_expVals_probeLevel_file = "/work/DAT_118__AML//Analysis//dset1/mRNA/raw_expVals_probeLevel.txt"
write.table(expVals,raw_expVals_probeLevel_file, row.names=T, col.names=T, quote=F)

#push the raw expression values to synapse
syn_raw_expVals_file = synStore(File(raw_expVals_probeLevel_file,parentId = 'syn2354332'), activity=act)

#standard error of measurement
se.expVals <- se.exprs(BSData)
head(se.expVals)

fData(BSData)
head(pData(BSData))

randIDs <- sample(featureNames(BSData),2000)
heatmap(expVals[randIDs,])

?addFeatureData