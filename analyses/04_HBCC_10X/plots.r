
# Making depth plots for HBCC_82041 10X

#Load packages
library(ggplot2)
library(dplyr)


# Read in regional file
# Read in files for ONT
b <- read.table("regional_cov_all.txt", header=T)

# Separate out rows
only_transcript_intron <- c('intron8transcript')
only_transcript_intron <- subset(b,REGION %in% only_transcript_intron)

exon9 <- c('exon9')
exon9 <- subset(b,REGION %in% exon9)

exon8 <- c('exon8')
exon8 <- subset(b,REGION %in% exon8)

intron8 <- c('intron8')
intron8 <- subset(b,REGION %in% intron8)

intron8_minus_transcript <- c('intron8_minus_transcript')
intron8_minus_transcript <- subset(b,REGION %in% intron8_minus_transcript)

cbPalette <- c("#F0E442","#E69F00", "#D55E00","#0072B2","#56B4E9", "#009E73", "#CC79A7")





# Create bar plot for only intron 8 transcript
plot <- ggplot(only_transcript_intron, aes(x=factor(REGION), y= NORMALIZEDDEPTH, fill=SAMPLE))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', na.rm=TRUE) + xlab("HBCC_82040") + 
  ylab("Normalized Coverage (mean depth/cell counts)") + ggtitle("40 bp Intron 8 Transcript Region 10X GG Coverage") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) + scale_fill_manual(values=cbPalette) + guides(fill=guide_legend(title="Cell Type"))
plot 
# Save the plot to the output file
ggsave("onlyintron8.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "onlyintron8.png", "\n"))


# Create bar plot for only intron 8 transcript
plot <- ggplot(exon9, aes(x=factor(REGION), y= NORMALIZEDDEPTH, fill=SAMPLE))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', na.rm=TRUE) + xlab("HBCC_82040") + 
  ylab("Normalized Coverage (mean depth/cell counts)") + ggtitle("Exon 9 10X GG Coverage") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) + scale_fill_manual(values=cbPalette) + guides(fill=guide_legend(title="Cell Type"))

# Save the plot to the output file
ggsave("exon9.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "exon9.png", "\n"))


# Create bar plot for only intron 8 transcript
plot <- ggplot(exon8, aes(x=factor(REGION), y= NORMALIZEDDEPTH, fill=SAMPLE))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', na.rm=TRUE) + xlab("HBCC_82040") + 
  ylab("Normalized Coverage (mean depth/cell counts)") + ggtitle("Exon 8 10X GG Coverage") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) + scale_fill_manual(values=cbPalette) + guides(fill=guide_legend(title="Cell Type"))

# Save the plot to the output file
ggsave("exon8.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "exon8.png", "\n"))


# Create bar plot for only intron 8 transcript
plot <- ggplot(intron8, aes(x=factor(REGION), y= NORMALIZEDDEPTH, fill=SAMPLE))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', na.rm=TRUE) + xlab("HBCC_82040") + 
  ylab("Normalized Coverage (mean depth/cell counts)") + ggtitle("Intron 8 10X GG Coverage") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) + scale_fill_manual(values=cbPalette) + guides(fill=guide_legend(title="Cell Type"))

# Save the plot to the output file
ggsave("intron8.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "intron8.png", "\n"))


# Create bar plot for only intron 8 transcript
plot <- ggplot(intron8_minus_transcript, aes(x=factor(REGION), y= NORMALIZEDDEPTH, fill=SAMPLE))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', na.rm=TRUE) + xlab("HBCC_82040") + 
  ylab("Normalized Coverage (mean depth/cell counts)") + ggtitle("Intron 8 minus 40bp Transcript Region 10X GG Coverage") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) + scale_fill_manual(values=cbPalette) + guides(fill=guide_legend(title="Cell Type"))


# Save the plot to the output file
ggsave("intron8_minus_transcript.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "intron8_minus_transcript.png", "\n"))




# Whole gene bar plot
b <- read.table("cov_all_whole.txt", header=T)
# The palette with grey:
cbPalette <- c("#F0E442","#E69F00", "#D55E00","#0072B2","#56B4E9", "#009E73", "#CC79A7")

plot <- ggplot(b, aes(x=factor(SAMPLE), y= NORMALIZEDDEPTH, fill=SAMPLE))
plot <- plot + geom_col(position = 'dodge', na.rm = TRUE) + xlab("HBCC_82040") + 
  ylab("Normalized Coverage (mean depth/cell counts)") + ggtitle("GBA1 10X GG Coverage") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 17))  + scale_fill_manual(values=cbPalette) +
  guides(fill=guide_legend(title="Cell Type"))
plot

# Save the plot to the output file
ggsave("GBA_whole.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "GBA_whole.png", "\n"))


# Whole gene bar plot
b <- read.table("cov_all_whole_GBAP1.txt", header=T)
# The palette with grey:
cbPalette <- c("#F0E442","#E69F00", "#D55E00","#0072B2","#56B4E9", "#009E73", "#CC79A7")

plot <- ggplot(b, aes(x=factor(SAMPLE), y= NORMALIZEDDEPTH, fill=SAMPLE))
plot <- plot + geom_col(position = 'dodge', na.rm = TRUE) + xlab("HBCC_82040") + 
  ylab("Normalized Coverage (mean depth/cell counts)") + ggtitle("GBAP1 10X GG Coverage") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 17)) +
scale_fill_manual(values=cbPalette) +
  guides(fill=guide_legend(title="Cell Type"))
plot

# Save the plot to the output file
ggsave("GBAP1_whole.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "GBAP1_whole.png", "\n"))
