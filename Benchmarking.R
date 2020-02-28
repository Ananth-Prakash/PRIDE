library(ggplot2)

## Compare iBAQ values between RUN1 (yoda; 10threads) and RUN2 (yoda; 33threads) PXD005819; (MQ v1.6.3.4)

RUN1 <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001608_10threads_noah/proteinGroups.txt', sep="\t", header=TRUE)
RUN2 <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001608_30threads_yoda/proteinGroups.txt', sep="\t", header=TRUE)

RUN1_iBAQ <- data.frame(RUN1[, c("Protein.IDs", "iBAQ")])
RUN1_iBAQ$Type <- rep("Linux_noah_10Threads", nrow(RUN1_iBAQ))

RUN2_iBAQ <- data.frame(RUN2[, c("Protein.IDs", "iBAQ")])
RUN2_iBAQ$Type <- rep("Linux_yoda_30Threads", nrow(RUN2_iBAQ))

Plotdata <- rbind(RUN1_iBAQ, RUN2_iBAQ)

ggplot(Plotdata, aes(x=iBAQ, colour = Type))+
  geom_density()+
  theme_bw()+
  scale_x_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD001608")

Merged_data <- merge(RUN1_iBAQ, RUN2_iBAQ,
                           by.x=c("Protein.IDs"), by.y=c("Protein.IDs"),
                           all.x=FALSE, all.y=FALSE)
#colnames(Merged_true_false_data)
colnames(Merged_data) <- c("Protein.IDs", "iBAQ.10Threads_noah", "Type.x", "iBAQ.30Threads_yoda", "Type.y")

ggplot(Merged_data, aes(x=iBAQ.10Threads_noah, y=iBAQ.30Threads_yoda))+
  geom_point(size = 1, alpha=0.2)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD001608")





### All benchmarking comparisions 
####################################

## Compare iBAQ values between Match_between_runs=False and Match_between_runs=True (MQ v1.6.3.4)

Match_bw_runs_false <- read.table('F:/PXD000666_Ananth_RefProt_MQv1.6.3.4_MatchBwRuns-False/combined/txt/proteinGroups.txt', sep="\t", header=TRUE)
Match_bw_runs_true <- read.table('F:/PXD000666_Ananth_RefProt_MQv1.6.3.4_MatchBwRuns-True/combined/txt/proteinGroups.txt', sep="\t", header=TRUE)

Match_bw_runs_false_iBAQ <- data.frame(Match_bw_runs_false[, c("Protein.IDs", "iBAQ")])
Match_bw_runs_false_iBAQ$Type <- rep("Match_between_runs_False", nrow(Match_bw_runs_false))

Match_bw_runs_true_iBAQ <- data.frame(Match_bw_runs_true[, c("Protein.IDs", "iBAQ")])
Match_bw_runs_true_iBAQ$Type <- rep("Match_between_runs_True", nrow(Match_bw_runs_true))

Plotdata <- rbind(Match_bw_runs_false_iBAQ, Match_bw_runs_true_iBAQ)

ggplot(Plotdata, aes(x=iBAQ, colour = Type))+
  geom_density()+
  theme_bw()+
  scale_x_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")


Merged_true_false_data <- merge(Match_bw_runs_false_iBAQ, Match_bw_runs_true_iBAQ,
                                by.x=c("Protein.IDs"), by.y=c("Protein.IDs"),
                                all.x=FALSE, all.y=FALSE)
#colnames(Merged_true_false_data)
colnames(Merged_true_false_data) <- c("Protein.IDs", "iBAQ.Match_bw_runs_False", "Type.x", "iBAQ.Match_bw_runs_True", "Type.y")

ggplot(Merged_true_false_data, aes(x=iBAQ.Match_bw_runs_False, y=iBAQ.Match_bw_runs_True))+
  geom_point(size = 1, alpha=0.2)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")



## Compare iBAQ values between Andy's output and Match_between_runs=True (MQ v1.6.3.4)

Andy_output <- read.table('F:/PXD000666_Andys_output/proteinGroups.txt', sep="\t", header=TRUE)
Andy_output_iBAQ <- data.frame(Andy_output[, c("Protein.IDs", "iBAQ")])
Andy_output_iBAQ$Type <- rep("Andy's_results", nrow(Andy_output_iBAQ))

Plotdata_ananth_andy <- rbind(Andy_output_iBAQ, Match_bw_runs_true_iBAQ)

ggplot(Plotdata_ananth_andy, aes(x=iBAQ, colour = Type))+
  geom_density()+
  theme_bw()+
  scale_x_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")


## Compare iBAQ values between Windows and Linux version PXD000666; Match between Runs =True (MQ v1.6.3.4)

Windows_VM_run <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse_PXD000666/RefProt_MQv1.6.3.4_MatchBwRuns-True/proteinGroups.txt', sep="\t", header=TRUE)
Linux_LSF_run <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse_PXD000666/RefProt_MQv1.6.3.4_Linux-LSF_MatchBwRuns-True/proteinGroups.txt', sep="\t", header=TRUE)

Windows_VM_iBAQ <- data.frame(Windows_VM_run[, c("Protein.IDs", "iBAQ")])
Windows_VM_iBAQ$Type <- rep("Windows_VM", nrow(Windows_VM_iBAQ))

Linux_LSF_iBAQ <- data.frame(Linux_LSF_run[, c("Protein.IDs", "iBAQ")])
Linux_LSF_iBAQ$Type <- rep("Linux_LSF", nrow(Linux_LSF_iBAQ))

Plotdata <- rbind(Windows_VM_iBAQ, Linux_LSF_iBAQ)

ggplot(Plotdata, aes(x=iBAQ, colour = Type))+
  geom_density()+
  theme_bw()+
  scale_x_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")


Merged_Windows_Linux_data <- merge(Windows_VM_iBAQ, Linux_LSF_iBAQ,
                                by.x=c("Protein.IDs"), by.y=c("Protein.IDs"),
                                all.x=FALSE, all.y=FALSE)
#colnames(Merged_true_false_data)
colnames(Merged_Windows_Linux_data) <- c("Protein.IDs", "iBAQ.Windows_VM", "Type.x", "iBAQ.Linux_LSF", "Type.y")

ggplot(Merged_Windows_Linux_data, aes(x=iBAQ.Windows_VM, y=iBAQ.Linux_LSF))+
  geom_point(size = 1, alpha=0.2)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")

## Compare iBAQ values between Linux (noah; 10threads) and Linux (version)yoda; 39threads) PXD000666; Match between Runs =True (MQ v1.6.3.4)

Linux_noah_10Threads_run <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/Mouse_PXD000666/RefProt_MQv1.6.3.4_Linux-noah_10Threads_MatchBwRuns-True_RUN1/proteinGroups.txt', sep="\t", header=TRUE)
Linux_yoda_39Threads_run <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/Mouse_PXD000666/RefProt_MQv1.6.3.4_Linux-yoda_39Threads_MatchBwRuns-True/proteinGroups.txt', sep="\t", header=TRUE)

Linux_noah_10Threads_iBAQ <- data.frame(Linux_noah_10Threads_run[, c("Protein.IDs", "iBAQ")])
Linux_noah_10Threads_iBAQ$Type <- rep("Linux_noah_10Threads", nrow(Linux_noah_10Threads_iBAQ))

Linux_yoda_39Threads_iBAQ <- data.frame(Linux_yoda_39Threads_run[, c("Protein.IDs", "iBAQ")])
Linux_yoda_39Threads_iBAQ$Type <- rep("Linux_yoda_39Threads", nrow(Linux_yoda_39Threads_iBAQ))

Plotdata <- rbind(Linux_noah_10Threads_iBAQ, Linux_yoda_39Threads_iBAQ)

ggplot(Plotdata, aes(x=iBAQ, colour = Type))+
  geom_density()+
  theme_bw()+
  scale_x_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")


Merged_Linux_data <- merge(Linux_noah_10Threads_iBAQ, Linux_yoda_39Threads_iBAQ,
                                   by.x=c("Protein.IDs"), by.y=c("Protein.IDs"),
                                   all.x=FALSE, all.y=FALSE)
#colnames(Merged_true_false_data)
colnames(Merged_Linux_data) <- c("Protein.IDs", "iBAQ.Linux_noah_10Threads", "Type.x", "iBAQ.Linux_yoda_39Threads", "Type.y")

ggplot(Merged_Linux_data, aes(x=iBAQ.Linux_noah_10Threads, y=iBAQ.Linux_yoda_39Threads))+
  geom_point(size = 1, alpha=0.2)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")


## Compare iBAQ values between Linux RUN1 (noah; 10threads) and Linux RUN2 (noah; 10threads) PXD000666; Match between Runs =True (MQ v1.6.3.4)

Linux_noah_RUN1_10Threads_run <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse_PXD000666/RefProt_MQv1.6.3.4_Linux-noah_10Threads_MatchBwRuns-True_RUN1/proteinGroups.txt', sep="\t", header=TRUE)
Linux_noah_RUN2_10Threads_run <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse_PXD000666/RefProt_MQv1.6.3.4_Linux-noah_10Threads_MatchBwRuns-True_RUN2/proteinGroups.txt', sep="\t", header=TRUE)

Linux_noah_RUN1_10Threads_iBAQ <- data.frame(Linux_noah_RUN1_10Threads_run[, c("Protein.IDs", "iBAQ")])
Linux_noah_RUN1_10Threads_iBAQ$Type <- rep("Linux_noah_10Threads_RUN1", nrow(Linux_noah_RUN1_10Threads_iBAQ))

Linux_noah_RUN2_10Threads_iBAQ <- data.frame(Linux_noah_RUN2_10Threads_run[, c("Protein.IDs", "iBAQ")])
Linux_noah_RUN2_10Threads_iBAQ$Type <- rep("Linux_noah_10Threads_RUN2", nrow(Linux_noah_RUN2_10Threads_iBAQ))

Plotdata <- rbind(Linux_noah_RUN1_10Threads_iBAQ, Linux_noah_RUN2_10Threads_iBAQ)

ggplot(Plotdata, aes(x=iBAQ, colour = Type))+
  geom_density()+
  theme_bw()+
  scale_x_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")


Merged_Linux_data <- merge(Linux_noah_RUN1_10Threads_iBAQ, Linux_noah_RUN2_10Threads_iBAQ,
                           by.x=c("Protein.IDs"), by.y=c("Protein.IDs"),
                           all.x=FALSE, all.y=FALSE)
#colnames(Merged_true_false_data)
colnames(Merged_Linux_data) <- c("Protein.IDs", "iBAQ.Linux_noah_10Threads_RUN1", "Type.x", "iBAQ.Linux_noah_10Threads_RUN2", "Type.y")

ggplot(Merged_Linux_data, aes(x=iBAQ.Linux_noah_10Threads_RUN1, y=iBAQ.Linux_noah_10Threads_RUN2))+
  geom_point(size = 1, alpha=0.2)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")

## Compare iBAQ values between Windows (RUN1) and Windows(RUN2) PXD000666; Match between Runs =True (MQ v1.6.3.4)

Windows_VM_run1 <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse_PXD000666/RefProt_MQv1.6.3.4_Windows-VM_MatchBwRuns-True_RUN1/proteinGroups.txt', sep="\t", header=TRUE)
Windows_VM_run2 <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse_PXD000666/RefProt_MQv1.6.3.4_Windows-VM_MatchBwRuns-True_RUN2/proteinGroups.txt', sep="\t", header=TRUE)

Windows_VM_run1_iBAQ <- data.frame(Windows_VM_run1[, c("Protein.IDs", "iBAQ")])
Windows_VM_run1_iBAQ$Type <- rep("Windows_VM_RUN1", nrow(Windows_VM_run1_iBAQ))

Windows_VM_run2_iBAQ <- data.frame(Windows_VM_run2[, c("Protein.IDs", "iBAQ")])
Windows_VM_run2_iBAQ$Type <- rep("Windows_VM_RUN2", nrow(Windows_VM_run2_iBAQ))

Plotdata <- rbind(Windows_VM_run1_iBAQ, Windows_VM_run2_iBAQ)

ggplot(Plotdata, aes(x=iBAQ, colour = Type))+
  geom_density()+
  theme_bw()+
  scale_x_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")


Merged_data <- merge(Windows_VM_run1_iBAQ, Windows_VM_run2_iBAQ,
                                   by.x=c("Protein.IDs"), by.y=c("Protein.IDs"),
                                   all.x=FALSE, all.y=FALSE)
#colnames(Merged_true_false_data)
colnames(Merged_data) <- c("Protein.IDs", "iBAQ.Windows_VM_RUN1", "Type.x", "iBAQ.Windows_VM_RUN2", "Type.y")

ggplot(Merged_data, aes(x=iBAQ.Windows_VM_RUN1, y=iBAQ.Windows_VM_RUN2))+
  geom_point(size = 1, alpha=0.2)+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  ggtitle("MaxQuant Benchmarking v1.6.3.4\nPXD000666")

