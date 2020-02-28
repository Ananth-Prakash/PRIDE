#To plot the differences in lengths of protein groups
library(ggplot2)

protein_group_length_diff <- read.table("/Users/ananth/Documents/PXD000666_MaxQuant_Benchmarking/RefProt_MQv1.6.3.4_MatchBwRuns-True/Protein_groups_length_difference_Andy_minus_RefProt_MQv1.6.3.4.txt",
                                        sep="\t", header=TRUE)

Median_length_diff <- median(protein_group_length_diff$Protein_groups_length_diff)
data_points <- nrow(protein_group_length_diff)
  
ggplot(protein_group_length_diff, aes(x=Protein_groups_length_diff))+
  geom_histogram(alpha=0.5, binwidth = 1, colour="black")+
  xlab("Protein groups length difference (Andy's result minus Match_bw_runs_True)")+
  geom_vline(xintercept = Median_length_diff, colour = "blue")+
  annotate("text", x = -25, y = 3000, label = paste("n =", data_points, sep = " "))+
  theme_bw()+
  ggtitle("Length difference bw protein groups\nMaxQuant Benchmarking v1.6.3.4 (PXD000666)")
