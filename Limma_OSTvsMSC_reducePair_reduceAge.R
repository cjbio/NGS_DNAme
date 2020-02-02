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
treatment_arrays <- factor(pData(norm_idat_filt)$Condition)
#batch_arrays <- factor(pData(norm_idat_filt)$Slide) #if differentt arraay
#batch_site <- factor(pData(norm_idat_filt)$Site) # if different collaborators
Sex <- factor(pData(norm_idat_filt)$predictedSex)
Age <- factor(pData(norm_idat_filt)$Sample_Group)
Pair <- factor(pData(norm_idat_filt)$pair)
design <- model.matrix(~0 + treatment_arrays + Pair)
colnames(design) <- colnames(design) %>% gsub("treatment_arrays", "",.)

treatment.design <- design[,1:2]
batch.design <- design[,-(1:2)]


## Batch corrected PCA ----

# remove batch effects for visualisation
batch_corr <- removeBatchEffect(getM(norm_idat_filt), 
                                covariates = batch.design,
                                design = treatment.design,
                                )
?removeBatchEffect()
pca <- batch_corr %>% t %>% prcomp



#save(norm_idat_filt, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/norm_idat_filt.Rdata")
#load("Projects/SRG/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis/Analysis2_results/norm_idat_filt.Rdata")
save(norm_idat_filt, file = "~/DNAme/norm_idat_filt.Rdata")

## Model fit ----

# M Value Selection
# Beta for fold change only

fitM <- lmFit(minfi::getM(norm_idat_filt), design)

contrasts <- makeContrasts(ostVsMsc = OST - MSC,
                           levels = design)
?c
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
  ttm <- topTable(fit2M, coef = i, number = Inf, p.value = pVal, adjust.method = "BH") %>% 
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

# quick view results
OSTvsMSC <- tt.out[[1]]
OSTvsMSC <- OSTvsMSC[order(OSTvsMSC$adj.P.Val),]
sum(OSTvsMSC$adj.P.Val < 0.05)
#sum(hipOAVsNOF$adj.P.Val < 0.05 & abs(hipOAVsNOF$Delta_Beta) > 0.1)
#sum(hipOAVsNOF$adj.P.Val < 0.05 & hipOAVsNOF$Delta_Beta > 0.1)
#sum(hipOAVsNOF$adj.P.Val < 0.05 & hipOAVsNOF$Delta_Beta < -0.1)
hist(OSTvsMSC$P.Value)



write.csv(OSTvsMSC[OSTvsMSC$adj.P.Val < 0.05,], file = "~/DNAme/sig_DMP_YoungvsOldOST_reduced_ageGroup.csv")
#            sep = "\t", row.names = FALSE, quote = FALSE)

#write.table(kneeOAVshipOA, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/kneeOAVshipOA_EPIC.txt",
#            sep = "\t", row.names = FALSE, quote = FALSE)


## Batch corrected beta values ----

#batch_corrBeta <- removeBatchEffect(getBeta(norm_idat_filt), 
#                                covariates = batch.design,
#                                design = treatment.design)

#rawBeta <- getBeta(norm_idat_filt)


toFind <- c("cg01154966", "cg20951255", "cg16474118", "cg00865429", "cg16683060", 
            "cg22021794", "cg14784820", "cg14770719", "cg18412185", "cg14359798", 
            "cg19447962", "cg01446731", "cg14063191", "cg11669284", "cg09709565", 
            "cg24249648", "cg08521995", "cg03656099")

#toFindBetas_batchCor <- batch_corrBeta[rownames(batch_corrBeta) %in% toFind,]
#colnames(toFindBetas_batchCor) <- as.character(sample_sheet$Sample_Name)

toFindBetas_rawBeta <- rawBeta[rownames(rawBeta) %in% toFind,]
colnames(toFindBetas_rawBeta) <- as.character(sample_sheet$Sample_Name)

write.table(rawBeta, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/rawBetas.txt",
            sep = "\t", quote = FALSE)

write.table(toFindBetas_rawBeta, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/rawBetas_subset18.txt",
            sep = "\t", quote = FALSE)

#write.table(batch_corrBeta, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/batchCorBetas.txt",
#            sep = "\t", quote = FALSE)

#write.table(toFindBetas_batchCor, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/batchCorBetas_subset18.txt",
#            sep = "\t", quote = FALSE)


## batch cor M values then convert to B values
## if batch correct directly on B values, negative values are generated

batchCorM2Beta <- m2beta(batch_corr)
colnames(batchCorM2Beta) <- as.character(targets$Sample.Names)
#toFindBetas_batchCorM2Beta <- batchCorM2Beta[rownames(batchCorM2Beta) %in% toFind,]
#colnames(toFindBetas_batchCorM2Beta) <- as.character(sample_sheet$Sample_Name)

write.table(batchCorM2Beta, file = "Projects/SRG/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis/Analysis_July2018/batchCorBetas.txt",
            sep = "\t", quote = FALSE)

#write.table(toFindBetas_batchCorM2Beta, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/batchCorBetas_subset18.txt",
#            sep = "\t", quote = FALSE)




########
write.table(annotation, file = "Documents/Work/Projects/SRG/DNA_methylation/Sulaco/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/probes_usedEPIC.txt",
            sep = "\t", quote = FALSE, row.names = F)


################################################################################
## calculate delta betas "by hand"
## July 2018

# load norm_idat_filt from code above

#pheno <- colData(norm_idat_filt)

#NOF_samples <- pheno[pheno$Type == "NOF",] %>% rownames(.)
#hip_samples <- pheno[pheno$Type == "hip",] %>% rownames(.)

#NOF_rawB <- rawBeta[,which(colnames(rawBeta) %in% NOF_samples)]
#NOF_rawB <- exp(rowMeans(log(NOF_rawB)))
#NOF_rawB <- NOF_rawB %>% rowMeans(.)

#hip_rawB <- rawBeta[,which(colnames(rawBeta) %in% hip_samples)]
#hip_rawB <- exp(rowMeans(log(hip_rawB)))
#hip_rawB <- hip_rawB %>% rowMeans(.)

#deltaB_NOF_hip <- hip_rawB - NOF_rawB

## Update fold change with batch corrected delta beta ----

batchB <- batchCorM2Beta

NOF_samples <- sample_sheet[sample_sheet$Type == "NOF",] %>% 
  .$Sample_Name %>% 
  droplevels(.)
hip_samples <- sample_sheet[sample_sheet$Type == "hip",] %>% 
  .$Sample_Name %>% 
  droplevels(.)
knee_samples <- sample_sheet[sample_sheet$Type == "knee",] %>% 
  .$Sample_Name %>% 
  droplevels(.)

NOF_batchB <- batchB[,which(colnames(batchB) %in% NOF_samples)] %>% rowMeans(.)
hip_batchB <- batchB[,which(colnames(batchB) %in% hip_samples)] %>% rowMeans(.)
knee_batchB <- batchB[,which(colnames(batchB) %in% knee_samples)] %>% rowMeans(.)

deltaB_NOF_hip <- hip_batchB - NOF_batchB %>% data.frame
deltaB_knee_hip <- knee_batchB - hip_batchB %>% data.frame

#hipOAVsNOF <- read.table("Projects/SRG/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/hipOAVsNOF_EPIC.txt",
#                         sep = "\t", header = TRUE, quote = "")

hipOAVsNOF <- hipOAVsNOF %>% merge(., deltaB_NOF_hip, 
                                   by.x = "CpG", 
                                   by.y = "row.names")

hipOAVsNOF$Delta_Beta <- hipOAVsNOF$.
hipOAVsNOF <- hipOAVsNOF[,-20]

write.table(hipOAVsNOF, file = "Projects/SRG/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis/Analysis_July2018/hipOAVsNOF_EPIC_batchBdelta.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

#kneeOAVshipOA <- read.table("Projects/SRG/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis2_results/kneeOAVshipOA_EPIC.txt",
#                            sep = "\t", header = TRUE, quote = "")

kneeOAVshipOA <- kneeOAVshipOA %>% merge(., deltaB_knee_hip, 
                                         by.x = "CpG", 
                                         by.y = "row.names")

kneeOAVshipOA$Delta_Beta <- kneeOAVshipOA$.
kneeOAVshipOA <- kneeOAVshipOA[,-20]

write.table(kneeOAVshipOA, file = "Projects/SRG/DNA_methylation/E181906_LouiseReynard_EPIC/Analysis/Analysis_July2018/kneeOAVshipOA_EPIC_batchBdelta.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
