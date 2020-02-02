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

#select MSCs only
select <- grep("MSC", targets$Sample.Names)
targets <- targets[select,]

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
d <- pca$x %>% as.data.frame %>% add_rownames("Sentrix_Tag") %>% 
  left_join({pData(norm_idat_filt) %>% 
      as.data.frame %>% 
      mutate(Sentrix_Tag = paste0(Slide,"_",Array))}) %>% 
  mutate_at(vars(c(85, 98)), funs(factor)) %>% 
  as.data.frame
pcv <- round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)

# change colour = Type (condition) or colour = Slide (array slide)
gg <- ggplot(d, aes(PC1,PC2)) +
  geom_point(aes(colour = Type), size = 2.5) + 
  theme_bw() +
  ggtitle("PCA of Normalised M Values") +
  theme(axis.title.x    = element_text(size=15),
        axis.title.y    = element_text(size=15)) +
  xlab(label = paste0("PC (", pcv[1], "%)")) +
  ylab(label = paste0("PC (", pcv[2], "%)")) +
  geom_text(aes(label = Sample_Name))
print(gg)

## Model design ----




# Define batches and create model
treatment_arrays <- factor(pData(norm_idat_filt)$Sample_Group)
#batch_arrays <- factor(pData(norm_idat_filt)$Slide) #if differentt arraay
#batch_site <- factor(pData(norm_idat_filt)$Site) # if different collaborators
Sex <- pData(norm_idat_filt)$predictedSex
#Age <- norm_idat_filt$Age

design <- model.matrix(~0 + treatment_arrays + Sex)
colnames(design) <- colnames(design) %>% gsub("treatment_arrays", "",.)

treatment.design <- design[,1:2]
batch.design <- design[,-(1:2)]


## using sva  to predict batch
mod <- model.matrix(~treatment_arrays + Sex)
mod0 <- model.matrix(~Sex)

# use norm counts and let SVA estimate number of batches. 
svobj <- sva(getM(norm_idat_filt), mod, mod0)


# add sv to sample sheet
targets <- data.frame(targets, 
                          SV1 = svobj$sv[,1],
                          SV2 = svobj$sv[,2])
SV1 <-  svobj$sv[,1]
SV2 <-  svobj$sv[,2]
# remove batch effect due to sample, FOR VISUALISATION ONLY
design <- model.matrix(~0 + treatment_arrays + Sex + SV1 + SV2)
colnames(design) <- colnames(design) %>% gsub("treatment_arrays", "",.)

treatment.design <- design[,1:2]
batch.design <- design[,-(1:2)]
## Batch corrected PCA ----

# remove batch effects for visualisation
batch_corr <- removeBatchEffect(getM(norm_idat_filt), 
                                covariates = batch.design,
                                design = treatment.design)

pca <- batch_corr %>% t %>% prcomp
d <- pca$x %>% as.data.frame %>% add_rownames("Sentrix_Tag") %>% 
  left_join({pData(norm_idat_filt) %>% 
      as.data.frame %>% 
      mutate(Sentrix_Tag = paste0(Slide,"_",Array))}) %>% 
  mutate_at(vars(c(39, 52)), funs(factor)) %>% 
  as.data.frame
pcv <- round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)

# change colour = Type (condition) or colour = Slide (array slide)
gg <- ggplot(d, aes(PC1,PC2)) +
  geom_point(aes(colour = Type), size = 2.5, alpha = 0.8) + 
  theme_bw() +
  ggtitle("PCA of Batch Corrected Normalised M Values") +
  theme(axis.title.x    = element_text(size=15),
        axis.title.y    = element_text(size=15)) +
  xlab(label = paste0("PC (", pcv[1], "%)")) +
  ylab(label = paste0("PC (", pcv[2], "%)")) +
  xlim(-250, 350) +
  geom_text(aes(label = Sample_Name)) +
  #scale_colour_manual(values = c("red", "blue"), 
  #                    labels = c("OA hip", "NOF")) +
  guides(color = guide_legend(title="Condition"))
print(gg)


#save(norm_idat_filt, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/norm_idat_filt.Rdata")
#load("Projects/SRG/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis/Analysis2_results/norm_idat_filt.Rdata")
save(norm_idat_filt, file = "~/DNAme/norm_idat_filt.Rdata")

## Model fit ----

# M Value Selection
# Beta for fold change only

fitM <- lmFit(minfi::getM(norm_idat_filt), design)

contrasts <- makeContrasts(youngVsOld = young - old,
                           levels = design)

#contrasts <- makeContrasts(hipOAVsNOF = hip - NOF,
#                           kneeOAVshipOA = knee - hip,
#                          levels = design)

fit2M <- eBayes(contrasts.fit(fitM, contrasts = contrasts))


# Get results
cont_names             <- colnames(contrasts)
pVal                  <- 1
deltab                <- 0

## Results ----

# select annotation columns
tt.out <- list()
for(i in cont_names) {
  print(i)
  ttm <- topTable(fit2M, coef = "treatment_arraysyoung", number = Inf, p.value = pVal, adjust.method = "BH") %>% 
    as.data.frame %>% add_rownames("CpG")
  if(nrow(ttm) > 0) {
    ttm  <- ttm %>% dplyr::select(CpG,P.Value,`adj.P.Val`)
    tt <- ttm %>% 
      left_join(annotation[,c(1:4, 18:19, 22:24, 32:37, 44:45)], by = c("CpG" = "Name")) 
    print(nrow(tt))
    tt.out[[i]]          <- tt
    #write_csv(tt, path = paste0("results/DM_",gsub(" - ","_-_",i),".csv"))
  } else {
    paste0(i, ": No DM Results") %>% print
  }
}
toptable(fit2M)


# quick view results
youngVsOld <- tt.out[[1]]
youngVsOld <- youngVsOld[order(youngVsOld$adj.P.Val),]
sum(youngVsOld$adj.P.Val < 0.05)
#sum(hipOAVsNOF$adj.P.Val < 0.05 & abs(hipOAVsNOF$Delta_Beta) > 0.1)
#sum(hipOAVsNOF$adj.P.Val < 0.05 & hipOAVsNOF$Delta_Beta > 0.1)
#sum(hipOAVsNOF$adj.P.Val < 0.05 & hipOAVsNOF$Delta_Beta < -0.1)
hist(youngVsOld$P.Value)



write.csv(youngVsOld[youngVsOld$adj.P.Val < 0.05,], file = "~/DNAme/sig_DMP_youngVsOld.csv")