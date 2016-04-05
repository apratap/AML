#source("http://bioconductor.org/biocLite.R")
#biocLite("crlmm")
#biocLite("humanomni5quadv1bCrlmm")
#source("http://bioconductor.org/biocLite.R")
#biocLite("VanillaICE")
#install.packages("doMC")

library("foreach")
library("crlmm")
library("humanomni5quadv1bCrlmm")
library(ff)
library(gdata)
library("VanillaICE")
library("parallel")
library(doMC)
registerDoMC(10)


#Ref: http://www.bioconductor.org/packages/2.13/bioc/vignettes/crlmm/inst/doc/IlluminaPreprocessCN.pdf

options(ffcaching="ffeachflush")
options("krlmm.cores" = 8)
outdir <- paste0("/work/DAT_118__AML//Analysis//dset1//CNV//analysis/", 'crlmm_output/')
ldPath(outdir)
dir.create(outdir)
ocProbesets(150e3)
ocSamples(500)

datadir <- "/work/DAT_118__AML/Analysis/dset1/CNV/data"

#creating a sample sheet as required by the package
# specific of sheet picked up from Appendix A
# ref : http://supportres.illumina.com/documents/myillumina/84861395-800d-4568-a38c-8efa6fda5015/genomestudio_gt_module_user_guide_reva.pdf
# reqd cols : Sample_ID, SentrixBarcode_A, SentrixPosition_A
#
master_list_file <- "/external-data/DAT_118__AML/dataset1-Raddich-lab/CARDINAL/metadata/Master_List.xls"
samplesheet <- read.xls(master_list_file,sheet=1,as.is=T, colClasses='character')
samplesheet <- samplesheet[,c('PatientID', 'SentrixPosition_A', 'Sentrix.Barcode', 'SNP.Array.Batch')]
colnames(samplesheet) <- c('Sample_ID','SentrixBarcode_A', 'SentrixPosition_A', 'Sample_Group')
#remove the rows with have ND
samplesheet <- samplesheet[ !samplesheet$SentrixBarcode_A == 'ND', ]

arrayNames <- paste(samplesheet$Sample_ID,samplesheet$SentrixPosition_A,samplesheet$SentrixBarcode_A,sep="_")

#sanity check
sum(file.exists(paste("/home/apratap/projects/AML/work/Analysis/dset1/CNV/data/",arrayNames, "_Grn.idat",sep="")))
sum(file.exists(paste("/home/apratap/projects/AML/work/Analysis/dset1/CNV/data/",arrayNames, "_Red.idat",sep="")))

#undebug(genotype.Illumina)
cnSet <- genotype.Illumina(sampleSheet = samplesheet,
                           arrayNames = arrayNames,
                           arrayInfoColNames = list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
                           path = datadir,
                           copynumber = T,
                           batch = samplesheet$Sample_Group,
                           cdfName = "humanomni5quadv1b",
                           verbose=T,
                           call.method='krlmm'
                          )
save(cnSet, file="~/projects//AML//work/Analysis/dset1//CNV//analysis/genotype_object.RData")
load(file="~/projects//AML//work/Analysis/dset1//CNV//analysis/genotype_object.RData")

#copynumber
cnSet.updated <- crlmmCopynumber(cnSet, type=c("SNP", "X.SNP") )
save(cnSet, file="~/projects//AML//work/Analysis/dset1//CNV//analysis/genotype_object.RData")

load(file="~/projects//AML//work/Analysis/dset1//CNV//analysis/genotype_object.RData")

library(VanillaICE)
open(cnSet)
oligoSetList <- BafLrrSetList(cnSet)
close(cnSet)
show(oligoSetList)
class(oligoSetList)
## oligoSnpSet of first chromosome
#oligoSetList[[1]]

#generate the log R ratio list
lrrList <- lrr(oligoSetList)
#generate the B allele frequency
bafList <- baf(oligoSetList)


library("plyr")
convert_to_df <- function(df){
  x <- as.data.frame(as.ffdf(df))
  x['probeID'] <- rownames(x) 
  x
}

#Log R ratios
logR_ratios <- ldply(.data=lrrList, .fun=convert_to_df)
rownames(logR_ratios) <- logR_ratios$probeID
logR_ratios$probeID <- NULL
colnames(logR_ratios) <- gsub('X', '', colnames(logR_ratios))
logR_ratios_file <- paste0(tempdir(), '/AML_dset1_CNV_logRratios.tsv')
write.table(logR_ratios, file=logR_ratios_file, sep="\t",row.names=T,quote=F )


# B Allele freq
Ballele_freq <- ldply(.data=bafList, .fun=convert_to_df)
rownames(Ballele_freq) <- Ballele_freq$probeID
Ballele_freq$probeID <- NULL
colnames(Ballele_freq) <- gsub('X', '', colnames(Ballele_freq))
head(Ballele_freq)
Ballele_freq_file <- paste0(tempdir(), '/AML_dset1_CNV_B_AlleleFreq.tsv')
write.table(Ballele_freq, file=Ballele_freq_file, sep="\t",row.names=T, quote=F )


#store the data on synapse
library("synapseClient")
synapseLogin()
synStore(File(logR_ratios_file, parentId='syn2785711'),
             used = 'syn2354339')




###########
##QC checks
###########
probeNames <- unlist(lapply(lrrList, rownames))
lapply(lrrList, colnames)

probeNames
length(probeNames)

calls(cnSet)[1:10,]

#SNR plot
library(lattice)
open(cnSet$SNR)

cnSet@phenoData

Ns(cnSet, i=1:3, j=3:4)

histogram(~cnSet$SNR[],breaks=100)

sampleNames(cnSet)
getHapMapIds(cnSet)

batch(cnSet)

ls(assayData(cnSet))

as.matrix(A(cnSet))[1:5,]

cnSet$SNR[]

hist
histogram(cnSet$SNR)
varLabels(cnSet)

print(histogram(~cnSet$SNR, panel=function(...){panel.histogram(...)},
                breaks=25, xlim=range(cnSet$SNR), xlab="SNR"))

# cnSet
# 
# ?genotype.Illumina
# 
# getCrlmmAnnotationName("humanomni5quadv1b")
# 
# samplesheet[10:12,]
# arrayNames[10:12]
# samplesheet$Sample_Group[10:12]
# 
# 
# samplesheet[20:30,]
# 
# genotype.Illumina
# 
# dim(samplesheet)
# arrayNames


load("~/projects//AML//work/Analysis/dset1//CNV//analysis/genotype_object.RData")

print(cnSet)