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
BiocManager::install("mygene")
library(mygene)


library(stringr)
library(dplyr)


setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012203_Ananth/')


#### 1. Read proteinGroups to get proteinID information.
proteinGroups_inp  <- read.table( "proteinGroups-SwissProt_only.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#### 2. Read evidence.txt file
evidence_file  <- read.table( "evidence-SwissProt_only.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#### 3. Read annotation file
#### Note: In the user defined the annotation file the RAW file names in coloumn 'Run' should not have .raw extension!
####       the Channel coloumn is case sensitive, and entries must only be in the format 'channel.1', 'channel.2', 'channel.3', ......
annot <- read.table("MSstat_annot.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")



#### I . Convert input files into MSstat format and do initial filtering
input.mq <- MaxQtoMSstatsTMTFormat(evidence_file, proteinGroups_inp, annot)

#### II Protein summarization.
quant.msstats <- proteinSummarization(input.mq,
                                      method="msstats")

#### III Test for all the possible pairs of conditions
#### For pairwise comparisons among all conditions: 
# test.pairwise <- groupComparisonTMT(quant.msstats)

#### If comparisons have to be made between certain conditions then specify
#### Check how many conditions are there to be compared. 
levels(quant.msstats$Condition)
# [1] "Alzheimers"   "Non_demented"

#### The numerator is set to 1 and denominator is set to -1
#### All other places in the matrix related to other conditions are set to 0
#### Here for example Alzheimers is set to denominator and hence -1
comparison<-matrix(c(-1,1),nrow=1)

#### Set the column names. (Has to be in the same order as levels(quant.msstats$Condition) )
colnames(comparison) <- levels(quant.msstats$Condition)

#### Set the row name for the matrix
row.names(comparison)<-"Non_demented / Alzheimers"



test.comparison <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison)


#### IV. Add gene name and protein name

#### 1. If the proteinGroups.txt file already has Protein.names and Gene.names annotated, use this section
#if( "Protein.names" %in% colnames(proteinGroups_inp) ){
# test.comparison.names <- merge(x = proteinGroups_inp[,c("Protein.IDs", "Protein.names", "Gene.names")], y = test.comparison,
#                               by.x=c("Protein.IDs"), by.y=c("Protein"),
#                               all.x=FALSE, all.y=TRUE)
#} else {
#### 2. Else add protein names and gene names
#### perform protein group to gene mapping
#### this loop will take some time

  temp <- test.comparison

  temp$"Gene.Symbol" <- "NA"
  temp$"Gene.Name" <- "NA"
  temp$"ENSG" <- "NA"
  temp$"unique.gene.count" <- "NA"

  for(i in 1:nrow(temp))
  {
    x <- data.frame(strsplit(as.character(temp[ i, "Protein"]), split = ";"), stringsAsFactors = FALSE)
  
    x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
    x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
    x[,1] <- gsub("-\\d+", "", x[,1], perl=TRUE) # remove isoform number
                  
    f = file()
    sink(file=f)
    res <- tryCatch(queryMany(x[,1], scopes="uniprot", fields=c("ensembl.gene", "symbol", "name"), species="human"), error = function(e) {print(0)})
    sink()
    close(f)

    print(res)
    
    if (class(res)=="DFrame")
    {
      temp[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
      temp[ i, "Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")
      temp[ i, "Gene.Name"] <- paste( unique(unlist(res$name[!is.na(res$name)])), collapse = ";") 
      #print(temp[i,"Gene.Symbol"])
      temp_symb <- temp[i,"Gene.Symbol"]
      temp[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
    }
  } #end of for
#} #end of else

test.comparison.names <- temp

# The output (test.comparison.names) has multiple entries of same proteins because of their association with multiple protein groups.
# To filter out duplicate entries and to retain only one protein group, we retain that entry which has all the proteins that are also the majority proteins
# ie., those proteins that have > 50% of the peptides belonging to that protein group.
deduped_proteingroups <- merge(x=test.comparison.names, y=data.frame(Majority.protein.IDs = proteinGroups_inp[,c('Majority.protein.IDs')]),
                               by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                               all.x=FALSE, all.y=FALSE) 
  
test.comparison.names.genecount.pvalue.filtered <- deduped_proteingroups[deduped_proteingroups$adj.pvalue <= 0.05 & deduped_proteingroups$unique.gene.count == 1,]

write.table(test.comparison.names, "Protein_log2FoldChange_pvalue.txt", sep = "\t", row.names = FALSE, quote = FALSE )
write.table(deduped_proteingroups, "Protein_log2FoldChange_pvalue-deduped_proteingroups.txt", sep = "\t", row.names = FALSE, quote = FALSE )
write.table(test.comparison.names.genecount.pvalue.filtered, "Protein_log2FoldChange_pvalue-filtered.txt", sep = "\t", row.names = FALSE, quote = FALSE )
