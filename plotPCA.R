
source("https://bioconductor.org/biocLite.R")
biocLite("methyAnalysis", )
biocLite("ChAMP", )
library(lumi)

colnames(beta) <- c("young2", "old7", "young1", "old3", "old4", "old5", "old6", "young3", "young8", "young7")
p16PCA <- na.omit(beta)
class(p16PCA)
pca= prcomp(p16PCA)
plot(pca$rotation[,1],pca$rotation[,2], xlab = "PC1", ylab = "PC2")
text(pca$rotation[,1],pca$rotation[,2], row.names(pca$rotation), cex=0.6, pos=1)



swan <- preprocessSWAN(RGSet, mSet = MSet, verbose = FALSE)
gset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                    removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                    quantileNormalize = TRUE, stratified = TRUE,
                                    mergeManifest = FALSE, sex = NULL)
gset.funnorm <- preprocessFunnorm(RGSet)


beta <- getBeta(ratioSet)
beta <- getBeta(swan)
beta <- getBeta(gset.quantile)
beta <- getBeta(gset.funnorm)
write.csv(beta, "betaFunnorm.csv")