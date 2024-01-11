## Making depth plots for HBCC_82041 10X
#Load packages
library(ggplot2)
library(dplyr)

cbPalette <- c("#56B4E9", "#CC79A7",  "#D55E00","#F0E442",  "#F0E442", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

# Read in files for ONT
b <- read.table("regional_cov_all.txt", header=T)

# Create new table with averages by genotype
GG <- b[b$GENOTYPE=='GG',]
GT <- b[b$GENOTYPE=='GT',]
TT <- b[b$GENOTYPE=='TT',]
GG_mean <- GG%>%
  group_by(REGION) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD=sd(NORMALIZEDDEPTH, na.rm = TRUE))
GG_mean <- GG_mean%>%
  mutate(GENOTYPE='GG')
GT_mean <- GT%>%
  group_by(REGION) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD= sd(NORMALIZEDDEPTH, na.rm = TRUE))
GT_mean <- GT_mean%>%
  mutate(GENOTYPE='GT')
TT_mean <- TT%>%
  group_by(REGION) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD= sd(NORMALIZEDDEPTH, na.rm = TRUE))
TT_mean <- TT_mean%>%
  mutate(GENOTYPE='TT')
genos_all <- rbind(GG_mean,GT_mean, TT_mean)


# Separate out rows



only_transcript_intron <- c('intron8transcript')
only_transcript_intron <- subset(genos_all,REGION %in% only_transcript_intron)

exon9 <- c('exon9')
exon9 <- subset(genos_all,REGION %in% exon9)

exon8 <- c('exon8')
exon8 <- subset(genos_all,REGION %in% exon8)

intron8 <- c('intron8')
intron8 <- subset(genos_all,REGION %in% intron8)

intron8_minus_transcript <- c('intron8_minus_transcript')
intron8_minus_transcript <- subset(genos_all,REGION %in% intron8_minus_transcript)




# Create bar plot for only intron 8 transcript
plot <- ggplot(only_transcript_intron, aes(x=factor(REGION), MEANDEPTH, fill=GENOTYPE)) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  geom_errorbar(aes(ymin = MEANDEPTH - SD, ymax = MEANDEPTH + SD),
                position = position_dodge(0.9), width = 0.25) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("40bp Intron 8, CRISPR LCL ONT RNAseq") +
   theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +  
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot 

# Save the plot to the output file
ggsave("onlyintron8.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "onlyintron8.png", "\n"))


# Create bar plot for only exon 9
plot <- ggplot(exon9, aes(x=factor(REGION), MEANDEPTH, fill=GENOTYPE)) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  geom_errorbar(aes(ymin = MEANDEPTH - SD, ymax = MEANDEPTH + SD),
                position = position_dodge(0.9), width = 0.25) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Exon 9, CRISPR LCL ONT RNASeq") +
   theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16))  +  
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot 

# Save the plot to the output file
ggsave("exon9.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "exon9.png", "\n"))


# Create bar plot for only exon 8
plot <- ggplot(exon8, aes(x=factor(REGION), MEANDEPTH, fill=GENOTYPE)) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  geom_errorbar(aes(ymin = MEANDEPTH - SD, ymax = MEANDEPTH + SD),
                position = position_dodge(0.9), width = 0.25) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Exon 8, CRISPR LCL ONT RNAseq") +
   theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +  
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot 

# Save the plot to the output file
ggsave("exon8.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "exon8.png", "\n"))


# Create bar plot for only intron 8 
plot <- ggplot(intron8, aes(x=factor(REGION), MEANDEPTH, fill=GENOTYPE)) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  geom_errorbar(aes(ymin = MEANDEPTH - SD, ymax = MEANDEPTH + SD),
                position = position_dodge(0.9), width = 0.25) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Intron 8, CRISPR LCL ONT RNAseq") +
   theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +  
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot 

# Save the plot to the output file
ggsave("intron8.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "intron8.png", "\n"))



# Create bar plot for only intron 8 minus transcript
plot <- ggplot(intron8_minus_transcript, aes(x=factor(REGION), MEANDEPTH, fill=GENOTYPE)) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  geom_errorbar(aes(ymin = MEANDEPTH - SD, ymax = MEANDEPTH + SD),
                position = position_dodge(0.9), width = 0.25) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Intron 8 minus 40bp, CRISPR LCL ONT RNAseq") +
   theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16))  +  
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot 

# Save the plot to the output file
ggsave("intron8_minus_transcript.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "intron8_minus_transcript.png", "\n"))




# Read in files for ONT whole gene
a <- read.table("cov_all_whole.txt", header=T)

# Create new table with averages by genotype
GG <- a[a$GENOTYPE=='GG',]
GT <- a[a$GENOTYPE=='GT',]
TT <- a[a$GENOTYPE=='TT',]
GG_mean <- GG%>%
  group_by(GENOTYPE) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD=sd(NORMALIZEDDEPTH, na.rm = TRUE))
GG_mean <- GG_mean%>%
  mutate(GENOTYPE='GG')
GT_mean <- GT%>%
  group_by(GENOTYPE) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD= sd(NORMALIZEDDEPTH, na.rm = TRUE))
GT_mean <- GT_mean%>%
  mutate(GENOTYPE='GT')
TT_mean <- TT%>%
  group_by(GENOTYPE) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD= sd(NORMALIZEDDEPTH, na.rm = TRUE))
TT_mean <- TT_mean%>%
  mutate(GENOTYPE='TT')
genos_all <- rbind(GG_mean,GT_mean, TT_mean)

# Create bar plot for averages across genomes
plot <- ggplot(genos_all, aes(GENOTYPE, MEANDEPTH, fill=GENOTYPE))
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + labs(x = NULL) + geom_errorbar(aes(ymin = MEANDEPTH - SD, ymax = MEANDEPTH + SD),
                position = position_dodge(0.9), width = 0.25) +
  ylab("Normalized Coverage") + ggtitle("GBA1, CRISPR LCL ONT RNAseq") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) + guides(fill=guide_legend(title="Genotype")) + scale_fill_manual(values = cbPalette)
plot


# Save the plot to the output file
ggsave("GBA_whole.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "GBA_whole.png", "\n"))

# Read in GBAP1 file
a <- read.table("cov_all_whole.GBAP1.txt", header=T)


# Create new table with averages by genotype
GG <- a[a$GENOTYPE=='GG',]
GT <- a[a$GENOTYPE=='GT',]
TT <- a[a$GENOTYPE=='TT',]
GG_mean <- GG%>%
  group_by(GENOTYPE) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD=sd(NORMALIZEDDEPTH, na.rm = TRUE))
GG_mean <- GG_mean%>%
  mutate(GENOTYPE='GG')
GT_mean <- GT%>%
  group_by(GENOTYPE) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD= sd(NORMALIZEDDEPTH, na.rm = TRUE))
GT_mean <- GT_mean%>%
  mutate(GENOTYPE='GT')
TT_mean <- TT%>%
  group_by(GENOTYPE) %>%
  summarise(MEANDEPTH=mean(NORMALIZEDDEPTH), SD= sd(NORMALIZEDDEPTH, na.rm = TRUE))
TT_mean <- TT_mean%>%
  mutate(GENOTYPE='TT')
genos_all <- rbind(GG_mean,GT_mean, TT_mean)


# Create bar plot for averages across genomes
plot <- ggplot(genos_all, aes(GENOTYPE, MEANDEPTH, fill=GENOTYPE))
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + geom_errorbar(aes(ymin = MEANDEPTH - SD, ymax = MEANDEPTH + SD),
                position = position_dodge(0.9), width = 0.25) + labs(x = NULL) + 
  ylab("Normalized Coverage") + ggtitle("GBAP1, CRISPR LCL ONT RNAseq") +
  theme_light() + theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) + guides(fill=guide_legend(title="Genotype")) + scale_fill_manual(values = cbPalette)
plot


# Save the plot to the output file
ggsave("GBAP1_whole.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "GBAP1_whole.png", "\n"))



# Read in files for ONT
b <- read.table("regional_cov_all.txt", header=T)

sample_order <- c("ND01137_GG_OG","ND01137_GG_Mock", "ND22789_GG", "ND01137_GT", "ND22789_GT", "ND22789_TT_Mock", "ND22789_TT_OG")
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





# Create bar plot for only intron 9 transcript
plot <- ggplot(only_transcript_intron, aes(x = factor(SAMPLE, levels = sample_order), y = NORMALIZEDDEPTH, fill = factor(GENOTYPE))) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("40bp Intron 8, CRISPR LCL ONT RNAseq") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.9, size=12),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +    # Center the title
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot


# Save the plot to the output file
ggsave("onlyintron8_persample.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "onlyintron8_persample.png", "\n"))


# Create bar plot for only exon9 transcript
plot <- ggplot(exon9, aes(x = factor(SAMPLE, levels = sample_order), y = NORMALIZEDDEPTH, fill = factor(GENOTYPE))) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Exon 9, CRISPR LCL ONT RNASeq") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.9, size=12),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +   # Center the title
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot

# Save the plot to the output file
ggsave("exon9_persample.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "exon9_persample.png", "\n"))



# Create bar plot for only exon 8 transcript
plot <- ggplot(exon8, aes(x = factor(SAMPLE, levels = sample_order), y = NORMALIZEDDEPTH, fill = factor(GENOTYPE))) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Exon 8, CRISPR LCL ONT RNAseq") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.9, size=12),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +    # Center the title
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot


# Save the plot to the output file
ggsave("exon8_persample.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "exon8_persample.png", "\n"))


# Create bar plot for only intron 8 transcript
plot <- ggplot(intron8, aes(x = factor(SAMPLE, levels = sample_order), y = NORMALIZEDDEPTH, fill = factor(GENOTYPE))) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Intron 8, CRISPR LCL ONT RNAseq") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.9, size=12),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +   # Center the title
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot


# Save the plot to the output file
ggsave("intron8_persample.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "intron8_persample.png", "\n"))


# Create bar plot for only intron 8 minus transcript
plot <- ggplot(intron8_minus_transcript, aes(x = factor(SAMPLE, levels = sample_order), y = NORMALIZEDDEPTH, fill = factor(GENOTYPE))) +
  geom_bar(stat = "identity", position = 'dodge', na.rm = TRUE) +
  labs(x = NULL) +
  ylab("Normalized Coverage") +
  ggtitle("Intron 8 minus 40bp, CRISPR LCL ONT RNAseq") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.9, size=12),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) +   # Center the title
  guides(fill = guide_legend(title = "Genotype")) +
  scale_fill_manual(values = cbPalette)

plot


# Save the plot to the output file
ggsave("intron8_minus_transcript_persample.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "intron8_minus_transcript_persample.png", "\n"))



# Read in files for ONT whole gene
a <- read.table("cov_all_whole.txt", header=T)


# Create bar plot for averages across genomes
plot <- ggplot(a, aes(x=factor(SAMPLE,levels=sample_order),y=NORMALIZEDDEPTH, fill = factor(GENOTYPE)))
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + labs(x = NULL) +
  ylab("Normalized Coverage") + ggtitle("GBA1, CRISPR LCL ONT RNAseq") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.9, size=12),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) + guides(fill=guide_legend(title="Genotype")) + scale_fill_manual(values = cbPalette)

plot


# Save the plot to the output file
ggsave("GBA_whole_persample.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "GBA_whole_persample.png", "\n"))

# Read in GBAP1 file
a <- read.table("cov_all_whole.GBAP1.txt", header=T)





# Create bar plot for averages across genomes
plot <- ggplot(a, aes(x=factor(SAMPLE,levels=sample_order),y=NORMALIZEDDEPTH, fill = factor(GENOTYPE)))
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + labs(x = NULL) + 
  ylab("Normalized Coverage") + ggtitle("GBAP1, CRISPR LCL ONT RNAseq") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.9, size=12),
        plot.title = element_text(hjust = 0.5, size=16), axis.title.y=element_text(size=16)) + guides(fill=guide_legend(title="Genotype")) + scale_fill_manual(values = cbPalette)
plot



# Save the plot to the output file
ggsave("GBAP1_whole_persample.png", plot, width = 7, height = 7)
cat(paste("Density plot saved as:", "GBAP1_whole_persample.png", "\n"))



