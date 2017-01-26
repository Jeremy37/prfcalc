#!/nfs/users/nfs_j/js29/R-3.0.3/bin/Rscript
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]

df <- read.delim(file)
epigenomeIDs <- colnames(df)[grep("_PRF", colnames(df))]
df.prf <- df[,epigenomeIDs]

df$MAX_PRF <- apply(df.prf, 1, function(x) max(x, na.rm=T))
df$Q90_PRF <- apply(df.prf, 1, function(x) quantile(x, probs=c(0.9), names=F, na.rm=T))
df$Q50_PRF <- apply(df.prf, 1, function(x) quantile(x, probs=c(0.5), names=F, na.rm=T))
df$Q10_PRF <- apply(df.prf, 1, function(x) quantile(x, probs=c(0.1), names=F, na.rm=T))
df$MIN_PRF <- apply(df.prf, 1, function(x) min(x, na.rm=T))
df$MAX_PRF[!is.finite(df$MAX_PRF)] = NA
df$MIN_PRF[!is.finite(df$MIN_PRF)] = NA

write.table(df, file="", quote=F, sep="\t", row.names=F)
