# This file has some common functions to calculate posterior probabilities
# of SNP association given a data.table object with P values and segment IDs.
# 
library(data.table)


getSegmentNaivePPAs = function(df)
{
  if (!("segment" %in% colnames(df))) {
    die("To get PPAs the data.table must have a 'segment' column.")
  }
  if (!("logBF" %in% colnames(df))) {
    if ("Pval" %in% colnames(df)) {
      # Get approximate Bayes Factors for each SNP
      df[,"logBF"] <- getLogBFFromPval(df)
    } else {
      die("To get PPAs the data.table must have either a Pval column or a logBF column.")
    }
  }
  segIndices <- getUniqueIndices(df$segment)
  ppas <- vector(mode="double")
  # For each gene/segment...
  for (i in 1:length(segIndices)) {
    startIndex <- segIndices[[i]][[2]]
    endIndex <- segIndices[[i]][[3]]
    ppas[startIndex:endIndex] <- getNaivePPA(df$logBF[startIndex:endIndex])
  }
  ppas
}

getNaivePPA = function(vecLogBF)
{
  logsegbfNaive <- -1000
  for (i in 1:length(vecLogBF)) {
    logsegbfNaive <- sumlog(logsegbfNaive, vecLogBF[i])
  }
  vecPPA <- vecLogBF - logsegbfNaive
  exp(vecPPA)
}

getLogBFFromPval = function(df)
{
  df[,'Z'] <- pToZStat(df$Pval)
  df$F[is.na(df$F)] <- 0.01
  df$N[is.na(df$N)] <- 1000
  apply(df[,c('Z','F','N'),with = FALSE], 1, function(d) calcLogBF(d['Z'],d['F'],d['N']))
}


getUniqueIndices = function(vec)
{
  if (class(vec) == "factor") {
    # For some reason factors are INCREDIBLY slow if used in the code below
    vec <- as.character(vec)
  }
  indicesList <- list()
  if (length(vec) == 0) {
    return(indicesList)
  }
  
  lastVal = vec[1]
  lastIndex = 1
  for (i in 1:length(vec)) {
    if (vec[i] != lastVal) {
      indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i-1))))
      lastVal <- vec[i]
      lastIndex <- i
    }
  }
  indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i))))
  return(indicesList)
}

pToZStat = function(pVals) {
  sqrt(qchisq(p=pVals, df=1, lower.tail=F))
}

calcLogBF = function(Z, f, N) {
  WW <- 0.1
  V <- approx_v(f, N)
  r <- WW / (V + WW)
  toreturn <- log( sqrt(1-r) ) + (Z*Z*r / 2)
  toreturn
}

approx_v = function(f, N) {
  1 / (2*f*(1-f) * N)
}

sumlog = function(logx, logy) {
  if (logx > logy) return(logx + log(1 + exp(logy-logx)))
  else return(logy + log(1 + exp(logx-logy)))
}


