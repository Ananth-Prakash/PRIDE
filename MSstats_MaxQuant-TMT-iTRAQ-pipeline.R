# MSstats tutorial for post-processing MaxQuant outputs - TMT/iTRAQ

# Webpage  - http://msstats.org/msstats-2/
# Tutorial - TMT/iTRAQ: https://www.bioconductor.org/packages/release/bioc/vignettes/MSstatsTMT/inst/doc/MSstatsTMT.html
# Tutorial - LabelFree: http://msstats.org/wp-content/uploads/2019/11/MSstats_v3.18.1_manual.pdf

# http://msstats.org/msstats-2/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstatsTMT")
library('MSstatsTMT', warn.conflicts = F, quietly = T, verbose = F)



# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#install.packages("remotes")
#remotes::install_github("Bioconductor/GenomeInfoDb")
#remotes::install_github("Bioconductor/GenomicFeatures")
BiocManager::install("mygene")
library(mygene)


library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD007160/')


#### 1. Read proteinGroups to get proteinID information.
proteinGroups_inp  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#### 2. Read evidence.txt file
evidence_file  <- read.table( "evidence-MS3.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#### 3. Read annotation file
#### Note: In the user defined the annotation file the RAW file names in coloumn 'Run' should not have .raw extension!
####       the Channel coloumn is case sensitive, and entries must only be in the format 'channel.1', 'channel.2', 'channel.3', ......
annot <- read.table("MSstat_annot_channels.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


##### 4. Read sample replicate information
##### Note: User defined annotation file
sample_replicates <- read.table("sample_replicates.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#####A. First we filter the proteinGroups input file before processing with MSstatsTMT
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "Reverse"] != "+", ]
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "Potential.contaminant"] != "+", ]
## consider those protein groups with more than 1 peptides matched
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "Peptides"] > 1, ]


# # Only consider those protein groups for which at least 50% of the sample replicates have non-zero intensity values

for(i in 1:nrow(sample_replicates))
{
  sample_replicates[i,] <- gsub("^","Reporter.intensity.corrected.", sample_replicates[i,], perl=TRUE)
  sample_replicates[i,] <- gsub(" ",".", sample_replicates[i,], perl=TRUE)
  
  subdata <- proteinGroups_inp[, c(unname(unlist(sample_replicates[i,-1])))]
  proteinGroups_inp$replicate_nonzero_instensity_count <- apply(subdata[1:ncol(subdata)], 1, function(x) length(which(x != 0)) )
  
  proteinGroups_inp <- proteinGroups_inp[proteinGroups_inp$replicate_nonzero_instensity_count >= ncol(subdata)/2,]
  
}


#####B. Add gene name and protein name
#### perform protein group to gene mapping
#### this loop will take some time

proteinGroups_inp$"ENSG" <- "NA"
proteinGroups_inp$"Gene.Symbol" <- "NA"
proteinGroups_inp$"Gene.Name" <- "NA"
proteinGroups_inp$"unique.gene.count" <- "NA"

for(i in 1:nrow(proteinGroups_inp))
{
  x <- data.frame(strsplit(as.character(proteinGroups_inp[ i, "Majority.protein.IDs"]), split = ";"), stringsAsFactors = FALSE)
  
  # x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
  # x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
  # x[,1] <- gsub("-\\d+", "", x[,1], perl=TRUE) # remove isoform number
  
  x[,1] <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
  
  f = file()
  sink(file=f)
  res <- tryCatch(queryMany(x[,1], scopes="uniprot", fields=c("ensembl.gene", "symbol", "name"), species="human"), error = function(e) {print(0)})
  sink()
  close(f)
  
  print(res)
  
  if (class(res)=="DFrame" | class(res) == "DataFrame")
  {
    proteinGroups_inp[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    proteinGroups_inp[ i, "Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")
    proteinGroups_inp[ i, "Gene.Name"] <- paste( unique(unlist(res$name[!is.na(res$name)])), collapse = ";") 
    #print(proteinGroups_inp[i,"Gene.Symbol"])
    temp_symb <- proteinGroups_inp[i,"Gene.Symbol"]
    proteinGroups_inp[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
  }
} #end of for


#####C Protein groups that are mapped to a unique gene
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp$ENSG != "" , ]
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp$ENSG != "NA" , ]
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "unique.gene.count"] == 1, ]
proteinGroups_inp <- proteinGroups_inp[ grep(";", proteinGroups_inp$ENSG, invert = TRUE) , ]

colnames(proteinGroups_inp)
# Although the first step of postprocessing  MaxQuant output in "process-protein-groups-with-gene-names.R"
# does aggregation, here it is not advised as it changes MaxQuant output format and MSstatsTMT does not recognise
# change in format.
# There will be slight change in number of rows before and after aggregation. 
# Normalised intensity values in examples such as TMT PXD007160 may slightly differ some entries.
## Aggregate the columns "Reporter.intensity.corrected". Check their column numbers
# data.to.map <- aggregate(data.to.map[ , 122:221], list("Gene ID" = data.to.map$ENSG, "Gene.Name" = data.to.map$Gene.Name), sum, na.rm =TRUE)



# Copy dataframe here, normalise in case of TMT using the "process-protein-groups-with-gene-names.R"
# script and use the normalised proteinGroups_inp as input to MSstatsTMT

#### Start MSstatsTMT
#### I . Convert input files into MSstat format and do initial filtering
input.mq <- MaxQtoMSstatsTMTFormat(evidence_file, proteinGroups_inp, annot,
                                   rmPSM_withMissing_withinRun = TRUE)


#### II Protein summarization.
quant.msstats <- proteinSummarization(input.mq,
                                      method="msstats",
                                      global_norm = TRUE,
                                      reference_norm = TRUE,
                                      remove_norm_channel = TRUE,
                                      remove_empty_channel = TRUE,)


#### III Test for all the possible pairs of conditions
#### For pairwise comparisons among all conditions: 
# test.pairwise <- groupComparisonTMT(quant.msstats)

#### If comparisons have to be made between certain conditions then specify
#### Check how many conditions are there to be compared. 
levels(quant.msstats$Condition)
# [1] "Alzheimer's disease","Alzheimer's disease/Parkinson's disease","Control","Global_internal_standard","Parkinson's disease" 

#### The numerator is set to 1 and denominator is set to -1
#### All other places in the matrix related to other conditions are set to 0

# comparison between "Alzheimer's disease" and "Control"
#### Here for example "Control" is set to denominator and hence -1
comparison_1 <- matrix(c(1,0,-1,0,0),nrow=1)
# comparison between "Alzheimer's disease/Parkinson's disease" and "Control"
comparison_2 <- matrix(c(0,1,-1,0,0),nrow=1)
# comparison between "Parkinson's disease" and "Control"
comparison_3 <- matrix(c(0,0,-1,0,1),nrow=1)

# comparison between "Alzheimer's disease" and "Alzheimer's disease/Parkinson's disease"
comparison_4 <- matrix(c(1,-1,0,0,0),nrow=1)
# comparison between "Alzheimer's disease" and "Parkinson's disease"
comparison_5 <- matrix(c(1,0,0,0,-1),nrow=1)
# comparison between "Parkinson's disease" and "Alzheimer's disease/Parkinson's disease"
comparison_6 <- matrix(c(0,-1,0,0,1),nrow=1)


#### Set the column names. (Has to be in the same order as levels(quant.msstats$Condition) )
colnames(comparison_1) <- levels(quant.msstats$Condition)
colnames(comparison_2) <- levels(quant.msstats$Condition)
colnames(comparison_3) <- levels(quant.msstats$Condition)
colnames(comparison_4) <- levels(quant.msstats$Condition)
colnames(comparison_5) <- levels(quant.msstats$Condition)
colnames(comparison_6) <- levels(quant.msstats$Condition)

#### Set the row name for the matrix
row.names(comparison_1)<-"AD / Control"
row.names(comparison_2)<-"AD-PD / Control"
row.names(comparison_3)<-"PD / Control"
row.names(comparison_4)<-"AD / AD-PD"
row.names(comparison_5)<-"AD / PD"
row.names(comparison_6)<-"PD / AD-PD"


test.comparison_1 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_1)
test.comparison_2 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_2)
test.comparison_3 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_3)
test.comparison_4 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_4)
test.comparison_5 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_5)
test.comparison_6 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_6)


# The output (test.comparison) has multiple entries of same proteins because of their association with multiple protein groups.
# To filter out duplicate entries and to retain only one protein group, we retain that entry which has all the proteins that are also the majority proteins
# ie., those proteins that have > 50% of the peptides belonging to that protein group.
deduped_proteingroups_1 <- merge(x=test.comparison_1, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                               by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                               all.x=FALSE, all.y=FALSE) 

deduped_proteingroups_2 <- merge(x=test.comparison_2, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                                 by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                                 all.x=FALSE, all.y=FALSE) 

deduped_proteingroups_3 <- merge(x=test.comparison_3, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                                 by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                                 all.x=FALSE, all.y=FALSE) 

deduped_proteingroups_4 <- merge(x=test.comparison_4, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                                 by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                                 all.x=FALSE, all.y=FALSE) 

deduped_proteingroups_5 <- merge(x=test.comparison_5, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                                 by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                                 all.x=FALSE, all.y=FALSE) 

deduped_proteingroups_6 <- merge(x=test.comparison_6, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                                 by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                                 all.x=FALSE, all.y=FALSE) 

group1 <- deduped_proteingroups_1[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group1) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "Alzheimer's Disease vs. Control.foldChange", 
                      "Alzheimer's Disease vs. Control.pValue")

group2 <- deduped_proteingroups_2[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group2) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "Co-morbid Alzheimer's & Parkinson's Disease vs. Control.foldChange",
                      "Co-morbid Alzheimer's & Parkinson's Disease vs. Control.pValue")

group3 <- deduped_proteingroups_3[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group3) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "Parkinson's Disease vs. Control.foldChange",
                      "Parkinson's Disease vs. Control.pValue")

group4 <- deduped_proteingroups_4[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group4) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "Alzheimer's Disease vs. Co-morbid Alzheimer's & Parkinson's Disease.foldChange",
                      "Alzheimer's Disease vs. Co-morbid Alzheimer's & Parkinson's Disease.pValue")

group5 <- deduped_proteingroups_5[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group5) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "Alzheimer's Disease vs. Parkinson's Disease.foldChange",
                      "Alzheimer's Disease vs. Parkinson's Disease.pValue")

group6 <- deduped_proteingroups_6[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group6) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "Parkinson's Disease vs. Co-morbid Alzheimer's & Parkinson's Disease.foldChange",
                      "Parkinson's Disease vs. Co-morbid Alzheimer's & Parkinson's Disease.pValue")

combined_foldchange_file <- merge(x=group1, y=group2,
                                  by.x=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  by.y=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  all.x=FALSE, all.y=FALSE) 

combined_foldchange_file <- merge(x=combined_foldchange_file, y=group3,
                                  by.x=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  by.y=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  all.x=FALSE, all.y=FALSE)

combined_foldchange_file <- merge(x=combined_foldchange_file, y=group4,
                                  by.x=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  by.y=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  all.x=FALSE, all.y=FALSE)

combined_foldchange_file <- merge(x=combined_foldchange_file, y=group5,
                                  by.x=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  by.y=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  all.x=FALSE, all.y=FALSE)

combined_foldchange_file <- merge(x=combined_foldchange_file, y=group6,
                                  by.x=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  by.y=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  all.x=FALSE, all.y=FALSE)

write.table(combined_foldchange_file, "FoldChange_groupcomparison.txt", sep = "\t", row.names = FALSE, quote = FALSE )


##### Additional processing
## This section looks at postprocessd output from 
# 1) FOT normalisation and 
# 2) MS-StatsTMT postprocessing after calculating log2 foldchange
# to keep same gene ids common between files.
# This setp is done as it is desired to submit to Expression Atlas in this format.

# NOTE: Postprocessing for FOT normalisation should be done and the output file be ready
# before runnng this step.

FOT_normalised_file <- read.table( "proteinGroups_final_GIS_BatchNormalised-ChannelsRenamed.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
log2_FC_file <- read.table( "FoldChange_groupcomparison.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
# Remove NA Pvalues
log2_FC_file <- log2_FC_file[complete.cases(log2_FC_file),]

FOT_genes <- FOT_normalised_file[,c("Gene.ID", "Gene.Symbol")]
log2_FC_genes <- log2_FC_file[,c("Gene.ID", "Gene.Name")]

# get common genes
merged_ids <- merge(x=FOT_genes, y=log2_FC_genes,
                    by.x=c("Gene.ID", "Gene.Symbol"), by.y=c("Gene.ID", "Gene.Name"),
                    all.x=FALSE, all.y=FALSE)

# filter these genes from both files
# the resulting file will have common genes between them.
FOT_normalised_file <- merge(x=merged_ids, y=FOT_normalised_file,
                             by.x=c("Gene.ID", "Gene.Symbol"),
                             by.y=c("Gene.ID", "Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)

log2_FC_file <- merge(x=merged_ids, y=log2_FC_file,
                             by.x=c("Gene.ID", "Gene.Symbol"),
                             by.y=c("Gene.ID", "Gene.Name"),
                             all.x=FALSE, all.y=FALSE)


write.table(FOT_normalised_file, "ExpressionAtlas_proteinGroups_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )
write.table(log2_FC_file, "ExpressionAtlas_FoldChange_groupcomparison.txt", sep = "\t", row.names = FALSE, quote = FALSE )
