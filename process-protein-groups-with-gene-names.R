## ---------------------------
##
## Script name: process-protein-groups.R 
##
## Purpose of script: process proteinGroups.txt file from MaxQuant to include in Expression Atlas
##
## Author: Dr. Andrew Jarnuczak
##
## Date Created: 2019-07-01
##
## Copyright (c) Andrew Jarnuczak, 2019
## Email: jarnuczak@ebi.ac.uk
##
## ---------------------------
##
## Notes: Steps performed in the script
## 1. clean up proteinGroups file (remove CONTAMINANTS and REVERSE)
## 2. for MS1 based quantification results (i.e. label free or SILAC) normalise iBAQ intensities to ppb
## 3. perform proteinGroup to Gene ID mapping based on mapping provided by UniProt 
## 4. perform any other necessary clean-up steps (this will depend on the dataset, for example if there are any "technical" assays included in the original dataset that should be removed before data is included in Expression Atlas) 
##   
##
## ---------------------------


# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")

library(mygene)

setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012203_Reanalysis/")

#In case that the library mygene was not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("mygene")

if ( !exists("EXP_TYPE")) warning("Please specify experiment type variable: EXP_TYE")
## Specify experiment type, i.e. what type of quantification is used, MS1-based or MS2-based. As a rule of thumb, label free and SILAC is MS1-based, and iTRAQ and TMT is MS2-based.
# EXP_TYPE <- "MS2-quant"
EXP_TYPE <- "MS2-quant"



##############
# define FOT normalisation function
# FOT stands for Fraction Of Total. In this normalisation method each protein iBAQ intensity value is scaled to the total amount of signal in a given MS run (column) and transformed to parts per billion (ppb)
fot.normalise <- function(x){
  data.sum <-   apply(x, 2, function(y){sum(y, na.rm=TRUE)})
  # barplot((data.sum), log = "y")
  #
  ##### do ppm normalisation
  x.mat <- as.matrix(x)
  x.mat.ppb <- apply(x.mat, 2, function(i) i/sum(i, na.rm = T) * 1000000000 )
  x.mat.ppb <- as.data.frame(x.mat.ppb)
  colnames(x.mat.ppb) <- paste("ppb.", colnames(x.mat.ppb), sep = "")
  return(x.mat.ppb)
}
##############


##### read the protein expression matrix produced by MaxQuant
tmp  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
## clean up 
tmp <- tmp[ tmp[, "Reverse"] != "+", ]
tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]
## consider those protein groups with more than 1 peptides matched
tmp <- tmp[ tmp[, "Peptides"] > 1, ]

## 
# if experiment is label free use iBAQ quant:
if(EXP_TYPE == "MS1-quant"){
  message("Collecting iBAQ quantification")
  tmp <- tmp[ ,c(2, grep("iBAQ.", colnames(tmp))) ]
}


# if experiment is TMT, use reporter intensities, if experiment iTRAQ use Intensities. Note, MaxQuant might change how it reports iTRAQ and TMT intensities.
if(EXP_TYPE == "MS2-quant"){
  message("Collecting MS2 intensities")
  if( any(grepl("Reporter.intensity.corrected", colnames(tmp))) ){
    
    # code modified by Ananth:
    # Filter to consider only protein groups for which the at least 50 % of the sample replicates 
    # have non-zero intensity values
    
    ##### Read sample replicate information
    ##### Note: User defined annotation file
    sample_replicates <- read.table("sample_replicates.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
    
    for(i in 1:nrow(sample_replicates))
    {
      sample_replicates[i,] <- gsub("^","Reporter.intensity.corrected.", sample_replicates[i,], perl=TRUE)
      
      subdata <- tmp[, c(unname(unlist(sample_replicates[i,-1])))]
      tmp$replicate_nonzero_instensity_count <- apply(subdata[1:ncol(subdata)], 1, function(x) length(which(x != 0)) )
      
      tmp <- tmp[tmp$replicate_nonzero_instensity_count >= ncol(subdata)/2,]
      
    }
    ### code modification stops
    tmp <- tmp[ ,c(2, grep("Reporter.intensity.corrected.[0-9].{1,}", colnames(tmp))) ]
  } else {
    tmp <- tmp[ ,c(2, grep("Intensity.", colnames(tmp))) ]
  }
}



#####
Majority.protein.IDs <- tmp$Majority.protein.IDs
tmp <- tmp[ , -1]
# for iTRAQ and TMT ppb normalization might not be the best method
if(EXP_TYPE == "MS1-quant"){
  tmp <- fot.normalise(tmp)  
}
#
tmp <- data.frame( cbind(Majority.protein.IDs, tmp, stringsAsFactors = FALSE) )
##############
tmp[tmp == 0] <- NA
tmp[ tmp == "NaN"] <- NA




  
##### perform the protein group to gene mapping
# this loop will take some time


data.to.map <- tmp
data.to.map$"ENSG" <- "NA"
data.to.map$"Gene.Name" <- "NA"

for(i in 1:nrow(data.to.map)){
    
    x <- data.frame(strsplit(data.to.map[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
    x_temp <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
    f = file()
    sink(file=f)
    a <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species="human"),
             error = function(e) {print(0)})
    sink()
    close(f)
                 # DFrame not DataFrame 
    if (class(a)=="DFrame"){
    data.to.map[ i, "ENSG"] <- paste( unique(unlist(a$ensembl.gene[!is.na(a$ensembl.gene)])), collapse = ";")
    data.to.map[ i, "Gene.Name"] <- paste( unique(unlist(a$symbol[!is.na(a$symbol)])), collapse = ";")}
    
    print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(data.to.map),1)),"%"))
  }



# remove protein groups that have no mapping to an ENSG gene IDs
data.to.map <- data.to.map[ data.to.map$ENSG != "" , ]
data.to.map <- data.to.map[ data.to.map$ENSG != "NA" , ]
#tmp <- data.to.map
# remove all protein groups that map to multiple ENSG gene IDs (this removes a lot of proteins) - the reasoning to remove these cases is that we cannot establish for sure which gene is contributing the signal to the protein abundance; all genes contribute equally or one gene is a majority? 
data.to.map <- data.to.map[ grep(";", data.to.map$ENSG, invert = TRUE) , ]
# for genes that map to multiple proteins, in order to determine the amount of protein that gene is producing we sum the protein quantification values
#xx.Majority.protein.IDs <- aggregate(data.to.map$Majority.protein.IDs, list(ESNG = data.to.map$ENSG ), function(x) paste0( (x) )  )
#
colnames(data.to.map)
## select which columns to aggregate
data.to.map <- aggregate(data.to.map[ , 2:(ncol(data.to.map)-2)], list("Gene ID" = data.to.map$ENSG, "Gene.Name" = data.to.map$Gene.Name), sum, na.rm =TRUE)
#
# data.to.map <- cbind(Majority.protein.IDs = (as.character(xx.Majority.protein.IDs$x)), data.to.map )


  
if(EXP_TYPE == "MS2-quant"){
  write.table(data.to.map, "proteinGroups_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )} else {
    write.table(data.to.map, "proteinGroups_ppb_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )   
  }

