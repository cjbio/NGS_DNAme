targets
## Packages ----
library(minfi)
library(DMRcate)
library(limma)
library(tidyverse)
#library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(pheatmap)
library(lumi)
rm(list=ls())
## Description ----
# Project: LR methylation data EPIC array- OA hip, OA knee and NOF
# Method: Analyse all samples together
# Method: Normalise using minfi functional normalisation only. No sva.
# Method: Remove SNPs, sex chromosomes, failed probes and cross reacting probes
# Date: analysis May 2018, updated beta fold changes July 2018

## UPDATE: Remove NOF sample RHH168 & , July 2018

## Prep pheno table ----

# Pheno data sample sheet, add Basename column for idats
targets <- read.metharray.sheet("~/DNAme/E192071_ChaoJiang_EPIC_220419")
targets$Basename<- paste(targets$Slide, targets$Array, sep = "_")

IDs <- read.csv("~/DNAme/DNA methylation seq samples (Chao).csv")

targets$Sample.Names <- IDs$Sample.Name
targets$Sample_Group <- IDs$AgeGroup
targets$Condition <- IDs$Condition

targets <- targets[c(1,5,7,12,8,13,9,14,10,15,11,4),]
targets$pair <-  c(1,1,2,2,3,3,4,4,5,5,6,6)
# read idats
RGset <- read.metharray.exp(targets = targets,
                            force = TRUE)

annotation <- getAnnotation(RGset)
annotation[annotation$Type == "I",] %>% rownames %>% length
annotation[annotation$Type == "II",] %>% rownames %>% length
(annotation[annotation$Type == "I",] %>% rownames %>% length) + 
  (annotation[annotation$Type == "II",] %>% rownames %>% length)
annoEPIC <- annotation

## QC probes ----
#filter step 1 - filter samples
lumi_dpval <- detectionP(RGset, type = "m+u") #detected pvalue of probe calls
lumi_passed <- lumi_dpval <= 0.05
call_rate <- colSums2(lumi_passed)/nrow(lumi_passed)
names(call_rate) <- colnames(lumi_dpval)

ggplot(data=data.frame(sample=sampleNames(RGset), call_rate=call_rate), 
       aes(x=sample, y=call_rate)) + 
  geom_point() +
  geom_hline(yintercept = 0.90, colour='red') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

lumi_failed <- !lumi_passed
sample_dpval_keep <- names(call_rate)[call_rate > 0.90] # remove samples with call_ratte < 0.9

## Normalise ----

# use functional normalisation method
raw_idat_filt <- RGset[,sample_dpval_keep]
norm_idat <- preprocessFunnorm(raw_idat_filt,
                               nPCs    = 2,
                               bgCorr  = TRUE,
                               dyeCorr = TRUE,
                               verbose = TRUE)

## Probe filters ----

# filter step  2 remove failed probes
Mvals <- minfi::getM(norm_idat)
remove <- which(Mvals==-Inf, arr.ind=TRUE)[,1]
if(length(remove) > 0) {
  norm_idat <- norm_idat[-remove,]
}else{
  norm_idat <- norm_idat
}
lumi_failed <- lumi_failed[,sample_dpval_keep]
probe_dpval_remove <- names(which(rowMeans(lumi_failed)>0.5, TRUE)) ##? 0.05
probe_dpval_remove_in <- probe_dpval_remove[!probe_dpval_remove %in% names(remove)]
norm_idat <- norm_idat[-match(probe_dpval_remove_in, rownames(norm_idat)),]

manifest <- getManifest(raw_idat_filt)

# remove SNPs (probes containning natural c/t polymorphism )
norm_idat_snp <- addSnpInfo(norm_idat)
norm_idat_snp_drop <- dropLociWithSnps(norm_idat_snp, 
                                       snps = c("SBE","CpG","Probe"), 
                                       maf  = 0.05)
annotation <- as.data.frame(getAnnotation(norm_idat_snp_drop))

# remove sex chromosomes
norm_idat_filt  <- norm_idat_snp_drop[-grep("chrX|chrY", annotation$chr),]

rm(Mvals, norm_idat_snp_drop, norm_idat_snp, manifest, remove)
gc()

# remove cross reacting probes (lists from LR)

CR1 <- read.csv("~/DNAme/13059_2016_1066_MOESM1_ESM.csv", header = TRUE)
#CR2 <- read.csv("Projects/SRG/DNA_methylation/crossreact-probes-Illumina450k (2).csv", header = TRUE)
CR_remove <- unique(as.character(CR1$CpG))

norm_idat_filt <- norm_idat_filt[-which(rownames(norm_idat_filt) %in% CR_remove),]

annotation <- as.data.frame(getAnnotation(norm_idat_filt))

## PCA ----

# mutate at group/type and slide variables. All samples - c(85, 98)
pca <- norm_idat_filt %>% getM %>% t %>% prcomp


## Model design ----

# Define batches and create model
treatment <- factor(pData(norm_idat_filt)$Condition)
#batch_arrays <- factor(pData(norm_idat_filt)$Slide) #if differentt arraay
#batch_site <- factor(pData(norm_idat_filt)$Site) # if different collaborators
#Sex <- factor(pData(norm_idat_filt)$predictedSex)
Age <- factor(pData(norm_idat_filt)$Sample_Group, levels = c("young", "old"))
#Pair <- factor(pData(norm_idat_filt)$pair)


design <- model.matrix(~Age+Age:treatment)

fit <- lmFit(minfi::getM(norm_idat_filt), design)
fit <- eBayes(fit)
colnames(fit)
topTable(fit, coef=3) # genes respond OST treatmeent in young donor
topTable(fit, coef=4) # genes respond OST treatmeent in old donor

fit2 <- contrasts.fit(fit, c(0,0,-1,1))
fit2 <- eBayes(fit2)
colnames(fit2)
topTable(fit2)
OSTvsMSC_age<- topTable(fit2, number = Inf)
OSTvsMSC_age <- OSTvsMSC_age[order(OSTvsMSC_age$adj.P.Val),]
sum(OSTvsMSC_age$adj.P.Val < 0.05)


################################################
design <- model.matrix(~Age*treatment)


fit <- lmFit(minfi::getM(norm_idat_filt), design)
cont.matrix <- cbind(SvsUinWT=c(0,0,1,0),SvsUinMu=c(0,0,1,1),Diff=c(0,0,0,1))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef = 1)

?contrasts.fit
sum(OSTvsMSC_age$adj.P.Val < 0.05)

topTable(fit2)
OSTvsMSC_age<- topTable(fit2, coef = 1, number = Inf)
OSTvsMSC_age <- OSTvsMSC_age[order(OSTvsMSC_age$adj.P.Val),]
sum(OSTvsMSC_age$adj.P.Val < 0.05)

##################################################
design <- model.matrix(~Age*treatment)
depfit <- duplicateCorrelation(minfi::getM(norm_idat_filt), design, block = Pair)

fit <- lmFit(minfi::getM(norm_idat_filt), design, correlation =depfit$consensus.correlation, block = Pair)
colnames(fit)


cont.matrix <- cbind(SvsUinWT=c(0,0,1,0),SvsUinMu=c(0,0,1,1),Diff=c(0,0,0,1))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef = 3)

?contrasts.fit
sum(OSTvsMSC_age$adj.P.Val < 0.05)

topTable(fit2)
OSTvsMSC_age<- topTable(fit2, coef = 3, number = Inf)
OSTvsMSC_age <- OSTvsMSC_age[order(OSTvsMSC_age$adj.P.Val),]
sum(OSTvsMSC_age$adj.P.Val < 0.05)

