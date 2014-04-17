library("crlmm")
library("humanomni5quadv1bCrlmm")
library(ff)
library(gdata)


#Ref: http://www.bioconductor.org/packages/2.13/bioc/vignettes/crlmm/inst/doc/IlluminaPreprocessCN.pdf

options(ffcaching="ffeachflush")
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


sum(file.exists(paste("/home/apratap/projects/AML/work/Analysis/dset1/CNV/data/",arrayNames, "_Grn.idat",sep="")))

sum(file.exists(paste("/home/apratap/projects/AML/work/Analysis/dset1/CNV/data/",arrayNames, "_Red.idat",sep="")))

undebug(genotype.Illumina)
cnSet <- genotype.Illumina(sampleSheet = samplesheet[1:5,],
                           arrayNames = arrayNames[1:5],
                           arrayInfoColNames = list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
                           path = datadir,
                           copynumber = T,
                           batch = samplesheet$Sample_Group[1:5],
                           cdfName = "humanomni5quadv1b",
                           call.method = "krlmm",
                           verbose=T
                          )


undebug(constructInf)
cnSet <- constructInf(sampleSheet = sampleSheet[1:5,], 
                      arrayNames = arrayNames[1:5],
                      path = datadir, 
                      arrayInfoColNames = arrayInfoColNames[1:5], 
                       XY = XY, 
                      cdfName = "humanomni5quadv1b",
                      verbose = T, 
                      batch =  samplesheet$Sample_Group[1:5], 
                      saveDate = T)

?crlmmIlluminaV2