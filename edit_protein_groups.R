# Script to remove sp|tr and accession annotations from protein group ids

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD005819_10threads/')

data <- read.table("proteinGroups.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
protein_ids <- data.frame(protein.ids = data[,1])

protein_ids$groups <- "NA"
for(i in 1:nrow(protein_ids)){
  x <- data.frame(strsplit(as.character(protein_ids[ i, "protein.ids"]), split = ";"), stringsAsFactors = FALSE)
  colnames(x) <- "prot"
  
  # Ananth: Edited to remove sp| or tr| and any characters that follow uniprot accessions.
  x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
  x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
  
  protein_ids[ i, "groups"] <- paste( x[,1], collapse = ";")
}

write.table(protein_ids, "proteinGroups_ids_only_without_sp_tr", sep = "\t", quote=FALSE, row.names = FALSE ) 
