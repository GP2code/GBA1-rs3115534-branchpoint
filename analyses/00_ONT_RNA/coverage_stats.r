a <- read.table("regional_cov_all.txt", header=T)

library(dplyr)

a <- a %>% 
  mutate(geno_new = case_when(
    GENOTYPE == "GG" ~ 1,
    GENOTYPE == "GT" ~ 1,
    GENOTYPE == "TT" ~ 0,
    TRUE ~ NA_integer_
  ))


# Separate out rows

only_transcript_intron <- c('intron8transcript')
only_transcript_intron <- subset(a,REGION %in% only_transcript_intron)

exon9 <- c('exon9')
exon9 <- subset(a,REGION %in% exon9)


exon8 <- c('exon8')
exon8 <- subset(a,REGION %in% exon8)

intron8 <- c('intron8')
intron8 <- subset(a,REGION %in% intron8)

intron8_minus_transcript <- c('intron8_minus_transcript')
intron8_minus_transcript <- subset(a,REGION %in% intron8_minus_transcript)

# Start glms
onlytranscript <- glm(NORMALIZEDDEPTH~geno_new, data=only_transcript_intron)
summary(onlytranscript)
model_summary <- summary(onlytranscript)
p_values <- model_summary$coefficients[, "Pr(>|t|)"]
print(p_values)

exon8glm <- glm(NORMALIZEDDEPTH~geno_new, data=exon8)
summary(exon8glm)
model_summary <- summary(exon8glm)
p_values <- model_summary$coefficients[, "Pr(>|t|)"]
print(p_values)


exon9glm <- glm(NORMALIZEDDEPTH~geno_new, data=exon9)
summary(exon9glm)
model_summary <- summary(exon9glm)
p_values <- model_summary$coefficients[, "Pr(>|t|)"]
print(p_values)

intron8glm <- glm(NORMALIZEDDEPTH~geno_new, data=intron8)
summary(intron8glm)
model_summary <- summary(intron8glm)
p_values <- model_summary$coefficients[, "Pr(>|t|)"]
print(p_values)

intronminusglm <- glm(NORMALIZEDDEPTH~geno_new, data=intron8_minus_transcript)
summary(intronminusglm)
model_summary <- summary(intronminusglm)
p_values <- model_summary$coefficients[, "Pr(>|t|)"]
print(p_values)



a <- read.table("cov_all_whole.txt", header=T)

a <- a %>% 
  mutate(geno_new = case_when(
    GENOTYPE == "GG" ~ 1,
    GENOTYPE == "GT" ~ 1,
    GENOTYPE == "TT" ~ 0,
    TRUE ~ NA_integer_
  ))

gba <- glm(NORMALIZEDDEPTH~geno_new, data=a)
summary(gba)
model_summary <- summary(gba)
p_values <- model_summary$coefficients[, "Pr(>|t|)"]
print(p_values)


a <- read.table("cov_all_whole.GBAP1.txt", header=T)

a <- a %>% 
  mutate(geno_new = case_when(
    GENOTYPE == "GG" ~ 1,
    GENOTYPE == "GT" ~ 1,
    GENOTYPE == "TT" ~ 0,
    TRUE ~ NA_integer_
  ))

gbap1 <- glm(NORMALIZEDDEPTH~geno_new, data=a)
summary(gbap1)
model_summary <- summary(gbap1)
p_values <- model_summary$coefficients[, "Pr(>|t|)"]
print(p_values)
