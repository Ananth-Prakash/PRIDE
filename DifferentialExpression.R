# Differential Expression Analysis for Proteomics
# https://bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DEP")

library("DEP")
library("dplyr")

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001608_10threads_yoda/')

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

proteinGroups  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
## clean up 
proteinGroups <- proteinGroups[ proteinGroups[, "Reverse"] != "+", ]
proteinGroups <- proteinGroups[ proteinGroups[, "Potential.contaminant"] != "+", ] 
## 
# if experiment is label free use iBAQ quant:
if(EXP_TYPE == "MS1-quant"){
  message("Collecting iBAQ quantification")
  proteinGroups <- proteinGroups[ ,c(2, grep("iBAQ.", colnames(proteinGroups))) ]
}


# if experiment is TMT, use reporter intensities, if experiment iTRAQ use Intensities. Note, MaxQuant might change how it reports iTRAQ and TMT intensities.
if(EXP_TYPE == "MS2-quant"){
  message("Collecting MS2 intensities")
  if( any(grepl("Reporter.intensity.corrected", colnames(proteinGroups))) ){
    proteinGroups <- proteinGroups[ ,c(2, grep("Reporter.intensity.corrected.[0-9].{1,}", colnames(proteinGroups))) ]
  } else {
    proteinGroups <- proteinGroups[ ,c(2, grep("Intensity.", colnames(proteinGroups))) ]
  }
}
#####
Majority.protein.IDs <- proteinGroups$Majority.protein.IDs
proteinGroups <- proteinGroups[ , -1]
# for iTRAQ and TMT ppb normalization might not be the best method
if(EXP_TYPE == "MS1-quant"){
  proteinGroups <- fot.normalise(proteinGroups)  
}
#
proteinGroups <- data.frame( cbind(Majority.protein.IDs, proteinGroups, stringsAsFactors = FALSE) )
##############
proteinGroups[proteinGroups == 0] <- NA
proteinGroups[ proteinGroups == "NaN"] <- NA


#######
# Perform UniProt protein ID to ensmbl gene IDs mapping
#######
# load the mapping reference file - this is a file provided by UniProt but can be any reference file 
# # download id mapping file for Human from UniProt
uniprot_url <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz"
#uniprot_url <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab.gz"
temp <- tempfile()
download.file(uniprot_url, temp, method = "libcurl", mode = "wb")
uniprot.map <- read.table(gzfile(temp), header = F, sep = "\t", fill = T, stringsAsFactors = FALSE)
unlink(temp)
colnames(uniprot.map) <- c("UniProtKB.AC","UniProtKB.ID","GeneID..EntrezGene.","RefSeq","GI", "PDB","GO","UniRef100","UniRef90","UniRef50","UniParc","PIR","NCBI.taxon","MIM","UniGene","PubMed","EMBL","EMBL.CDS","Ensembl","Ensembl_TRS","Ensembl_PRO","Additional.PubMed")


#################
# Add Gene names
#################

uniprot_url2 <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
temp <- tempfile()
download.file(uniprot_url2, temp, method = "libcurl", mode = "wb")
uniprot.map_gene_names <- read.table(gzfile(temp), header = F, sep = "\t", fill = T, stringsAsFactors = FALSE)
unlink(temp)
uniprot.map_gene_names <- uniprot.map_gene_names[uniprot.map_gene_names$V2 == "Gene_Name",]
colnames(uniprot.map_gene_names) <- c("UniProt_Acc", "Type", "Gene.name")

uniprot.map <- merge(uniprot.map, uniprot.map_gene_names,
                     by.x=c("UniProtKB.AC"), by.y=c("UniProt_Acc"),
                     all.x=TRUE, all.y=FALSE)

### Some UniProt Accessions have more than one gene names
# uniprot.map_gene_names[uniprot.map_gene_names$UniProt_Acc %>% duplicated == TRUE,]

uniprot.map <- uniprot.map[ , c(1, 19, 24)]



# a note about Accession and ID numbers in Uniprot: https://www.uniprot.org/help/difference_accession_entryname
# What is the difference between an accession number (AC) and the entry name?
#   
#   An accession number (AC) is assigned to each sequence upon inclusion into UniProtKB. Accession numbers are stable from release to release. If several UniProtKB entries are merged into one, for reasons of minimizing redundancy, the accession numbers of all relevant entries are kept. Each entry has one primary AC and optional secondary ACs.
# 
# The 'Entry name' (formerly ID) is a unique identifier, often containing biologically relevant information. It is sometimes necessary, for reasons of consistency, to change IDs (for example to ensure that related entries have similar names). Another common cause for changing an ID is when an entry is promoted from UniProtKB's TrEMBL section (with computationally-annotated records) to the Swiss-Prot section (with fully curated records). However, an accession number (AC) is always conserved, and therefore allows unambiguous citation of UniProt entries. If a UniProtKB entry contains more than one accession number, the first one (or primary accession number) should be cited. 




##### perform the protein group to gene mapping
# this loop will take some time
# we might want to rewrite this step to speed things up, maybe use API calls instead
# ?
data.to.map <- proteinGroups
data.to.map$ENSG <- "NA"
for(i in 1:nrow(data.to.map)){
  x <- data.frame(strsplit(data.to.map[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
  colnames(x) <- "prot"
  
  # Ananth: Edited to remove sp| or tr| and any characters that follow uniprot accessions.
  x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
  x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
  
  # extract canonical UniProt protein ID
  x[,1] <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
  x <- merge(x, uniprot.map, by.x = "prot", by.y = "UniProtKB.AC" )
  all.genes <- sub(" ","", unlist(strsplit(x[ ,2], ";" )))
  data.to.map[ i, "ENSG"] <- paste( unique(all.genes), collapse = ";")
}


# remove protein groups that have no mapping to an ENSG gene IDs
data.to.map <- data.to.map[ data.to.map$ENSG != "" , ]
# remove all protein groups that map to multiple ENSG gene IDs (this removes a lot of proteins) - the reasoning to remove these cases is that we cannot establish for sure which gene is contributing the signal to the protein abundance; all genes contribute equally or one gene is a majority? 
data.to.map <- data.to.map[ grep(";", data.to.map$ENSG, invert = TRUE) , ]
# for genes that map to multiple proteins, in order to determine the amount of protein that gene is producing we sum the protein quantification values
xx.Majority.protein.IDs <- aggregate(data.to.map$Majority.protein.IDs, list(ESNG = data.to.map$ENSG ), function(x) paste0( (x) )  )


data.to.map <- merge(data.to.map, uniprot.map,
                     by.x=c("ENSG"), by.y=c("Ensembl"),
                     all.x=TRUE, all.y=FALSE)


colnames(data.to.map)
## select which columns to aggregate

colnames(data.to.map)[14] <- "ENSG"
data_unique <- make_unique(data.to.map, "ENSG", "Majority.protein.IDs", delim = ";")

Reporter_intensity_columns <- grep("Reporter.intensity.corrected.", colnames(data.to.map)) # get Reporter.intensity.corrected column numbers
data_se_parsed <- make_se_parse(data.to.map, Reporter_intensity_columns)
