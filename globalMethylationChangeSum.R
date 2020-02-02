filtM <- getM(norm_idat_filt)
colnames(filtM) <- targets$Sample_Group

filtM <- as.data.frame(filtM)
hist(filtM$young)
hist(filtM$old)

test2 <- as.data.frame(
  sapply(unique(names(filtM)), 
         function(col) rowMeans(filtM[names(filtM) == col]) # calculate row means
  )
)
filtM[1,]
test2[1,]
4.169472+4.12602+3.712828+4.038363+4.193516
20.2402/5

hist(test2$young)
hist(test2$old)
sum(test2$young)
sum(test2$old)


mean(test2$young)
mean(test2$old)

table(test2)

markers <- read.csv("~/DNAme/runx2_motif_genes.csv", header = TRUE)

?getAnnotation()
View(youngVsOld)

fitM[rownames(fitM)=="cp09780241",]
