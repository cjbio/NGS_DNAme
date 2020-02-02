pca <- norm_idat_filt %>% getM %>% t %>% prcomp
pcv <- round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)
d <- as.data.frame(pca$x)
rownames(d) <- targets$Sample.Names


AgeGroup <- factor(targets$Sample_Group, levels = c("young", "old"))
colourByGroup <- sub("old", "black", targets$Sample_Group)
colourByGroup <- sub("young", "red", colourByGroup)


gg <- ggplot(d, aes(x = PC1,y = PC2, colour = AgeGroup)) +
  geom_point(size=3) +
  scale_colour_manual(values=colourByGroup) +
  #geom_point(color = "pink") +
  #geom_text(aes(label = rownames(d))) +
  xlab(label = paste0("PC1 (", pcv[1], "%)")) +
  ylab(label = paste0("PC2 (", pcv[2], "%)"))
  #geom_text(aes(label = pData(norm_idat_filt)$predictedSex)) #check for gender bias
print(gg)
round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)
#plot(pca$rotation[,1],pca$rotation[,2], xlab = "PC1", ylab = "PC2")
#text(pca$rotation[,1],pca$rotation[,2], row.names(pca$rotation), cex=0.6, pos=1)

############

pca= prcomp(norm_idat_filt %>% getBeta)
gg <- ggplot(as.data.frame(pca$rotation), aes(x = PC1,y = PC2, colour = AgeGroup)) +
  geom_point(size=3) +
  scale_colour_manual(values=colourByGroup) +
  #geom_point(color = "pink") +
  #geom_text(aes(label = rownames(d))) +
  xlab(label = paste0("PC1 (", pcv[1], "%)")) +
  ylab(label = paste0("PC2 (", pcv[2], "%)"))
#geom_text(aes(label = pData(norm_idat_filt)$predictedSex)) #check for gender bias
print(gg)
plot(pca$rotation[,1],pca$rotation[,2], xlab = "PC1", ylab = "PC2")
text(pca$rotation[,1],pca$rotation[,2], row.names(pca$rotation), cex=0.6, pos=1)
###########
library("calibrate")
## I like mine better:
rld_pca <- function (batch_corr, ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA plot", textcx=0.55, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  pca = batch_corr %>% t %>% prcomp
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), pch=20, cex = 2,col=colourByGroup, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  #with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  #legend(legendpos, legend=c("young","old"), col=c(20,18), pch=c(20))
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
#png("qc-pca-allpassages.png", 16000, 16000, pointsize=40)
rld_pca(batch_corr, xlim(-250, 350))
