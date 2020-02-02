library(wateRmelon)



betaValues <- getBeta(norm_idat_filt)
colnames(betaValues) <- targets$Sample.Names

length(rownames(betaValues))

ageMSCs <- agep(betaValues)

write.csv(ageMSCs, "HovMSCsPredictedAge.csv")
