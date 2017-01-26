#!/software/R-3.2.2/bin/Rscript
library(data.table)
library(RobustRankAggreg)
#source("ppaFunctions.R")

myargs <- NULL
epigenomeIDs <- NULL



main = function()
{
  myargs <<- getArgs()
  print(myargs)
  
  if (is.null(myargs$file)) {
    stop("Missing parameter --file")
  }
  if (is.null(myargs$verbose)) {
    myargs$verbose <<- F
  }
  
  if (is.null(myargs$metadata)) {
    myargs$metadata = "/lustre/scratch110/sanger/dg13/share/jeremy/prfscores/epigenomes/epigenome.metadata.txt"
  }
  metadata.df <- read.delim(myargs$metadata, stringsAsFactors=F, header=T)
  metadata.df <- metadata.df[order(metadata.df$Epigenome.ID),]
  rownames(metadata.df) <- metadata.df$Epigenome.ID
  
  options(datatable.fread.datatable=T)
  options(stringsAsFactors = T)
  options(datatable.verbose = F)
  options(datatable.showProgress = F)
  
  df <- fread(myargs$file, stringsAsFactors=T, header=T, sep='\t', na.strings=c("NA","N/A",""))
  
  if (!is.null(myargs$segmentCol)) {
    df$segment <- df[,myargs$segmentCol,with=F][[1]]
  }
  if (!is.null(myargs$ppaCol)) {
    df$PPANaive <- as.double(df[,myargs$ppaCol,with=F][[1]])
  }
  if (!is.null(myargs$pCol)) {
    df$Pval <- as.double(df[,myargs$pCol,with=F][[1]])
  }
  if (!is.null(myargs$F)) {
    df$F <- as.double(df[,myargs$F,with=F][[1]])
  }
  if (!is.null(myargs$NCol)) {
    df$N = as.integer(df[,myargs$NCol,with=F][[1]])
  } else if (!is.null(myargs$N)) {
    df$N = as.integer(myargs$N)
  }
  if (!is.null(myargs$debug)) {
    myargs$debug = T
  } else {
    myargs$debug = F
  }
  if (!is.null(myargs$perlocus)) {
    myargs$perlocus = T
  } else {
    myargs$perlocus = F
  }
  
  # Check that required columns are present
  if (! ("PPANaive" %in% colnames(df) | ("Pval" %in% colnames(df) & "N" %in% colnames(df) & "F" %in% colnames(df)))) {
    stop("Missing required column: either 'PPANaive' column or 'Pval', 'N', and 'F' columns must be present.")
  }
  
  if (is.null(myargs$epigenomes)) {
    write("Inferring epigenomes from header...", stderr())
    epigenomeIndexes = grep("_PRF", colnames(df))
    colnames(df)[epigenomeIndexes] <- gsub("_PRF","",colnames(df)[epigenomeIndexes])
    epigenomeIDs <- colnames(df)[epigenomeIndexes]
    write(paste0("Epigenome IDs: ", paste(epigenomeIDs, collapse=",")), stderr())
  } else {
    epigenomeIDs <- strsplit(myargs$epigenomes, ",", fixed=T)
    if (! all(epigenomeIDs %in% colnames(df))) {
      stop(paste0("Not all epigenome IDs in list (", myargs$epigenomes, ") are in columns of file ", myargs$file))
    }
  }
  
  if ("Pval" %in% colnames(df)) {
    # Get approximate Bayes Factors for each SNP
    df[,"logBF"] <- getBFFromPval(df)
  }
  if (!("PPANaive" %in% colnames(df))) {
    # Get the posterior probability of association for each SNP
    # Requires logBF to have been defined
    df[,"PPANaive"] <- getSegmentNaivePPAs(df)
  }
  
  
  if (myargs$debug) {
    write.table(df, "output.temp.txt", row.names=F, quote=F, sep="\t")
  }
  
  # Get weighted PRF scores for each segment for each epigenome
  segNames <- unique(df$segment)
  segIndices <- getUniqueIndices(df$segment)
  if (myargs$debug) {
    print(segNames)
    print(segIndices)
  } 
  weightedPRFMat <- matrix(nrow=length(epigenomeIDs), ncol=length(segIndices), dimnames=list(epigenomeIDs, segNames))
  weightedMaxPRFMat <- matrix(nrow=length(epigenomeIDs), ncol=length(segIndices), dimnames=list(epigenomeIDs, segNames))
  funcWeightedPRFMat <- matrix(nrow=length(epigenomeIDs), ncol=length(segIndices), dimnames=list(epigenomeIDs, segNames))
  maxPRFMat <- matrix(nrow=length(epigenomeIDs), ncol=length(segIndices), dimnames=list(epigenomeIDs, segNames))
  maxPPAMat <- matrix(nrow=length(epigenomeIDs), ncol=length(segIndices), dimnames=list(epigenomeIDs, segNames))
  
  for (i in 1:length(segIndices)) {
    startIndex <- segIndices[[i]][[2]]
    endIndex <- segIndices[[i]][[3]]
    # Get SNPs in this segment, and consider only those in the 95% statistical credible set
    df.seg <- df[startIndex:endIndex,]
    df.seg.cred <- df.seg[PPANaive > 0.001,][order(PPANaive, decreasing=T),]
    df.seg.cred <- df.seg.cred[, cumPPANaive := cumsum(PPANaive)]
    credSetInd <- which(df.seg.cred$cumPPANaive <= 0.95)
    if (length(credSetInd) == 0) {
      credSetMax <- 1
    } else {
      credSetMax <- max(credSetInd)
    }
    df.seg.cred <- df.seg.cred[1:credSetMax]
    
    # Get a "baseline" PRF score that respresents some normal PRF score
    # across epigenomes for this region. This is used to scale 
    prfBaseline = median(unlist(df.seg.cred[,epigenomeIDs,with=F]))
    #df.seg.cred[,epigenomeIDs,with=F] <- df.seg.cred[,epigenomeIDs,with=F] - prfBaseline
    for (j in 1:length(epigenomeIDs)) {
      epigenomeID <- epigenomeIDs[j]
      df.eid <- df.seg.cred[,c("PPANaive", epigenomeID, "logBF"),with=F]
      prfScores <- df.eid[,epigenomeID,with=F][[1]]

      if (all(is.na(prfScores))) {
        weightedPRFMat[j,i] <- NA
        weightedMaxPRFMat[j,i] <- NA
        maxPRFMat[j,i] <- NA
        funcWeightedPRFMat[j,i] <- NA
      } else {
        weightedPRFMat[j,i] <- sum(df.eid$PPANaive * prfScores, na.rm=T)
        weightedMaxPRFMat[j,i] <- max(df.eid$PPANaive * prfScores, na.rm=T)
        maxPRFMat[j,i] <- max(prfScores, na.rm=T)
        
        funcPPAs <- getFunctionalPPA(prfScores, df.eid$logBF)
        funcWeightedPRFMat[j,i] <- sum(funcPPAs * prfScores, na.rm=T)
      }
      
      prf <- df[startIndex:endIndex,epigenomeID,with=F][[1]]
      logBF <- df$logBF[startIndex:endIndex]
      # If there is only one element in the vectors, we add a dummy element
      # with ln(7) lower logBF. This prevents the PPA from being exactly 1
      # and so enables ranking differences in PPA based on the PRF for these
      # cases.
      if (length(logBF) == 1) {
        prf[2] <- 0
        logBF[2] <- logBF[1] - 7
      }
      maxPPAMat[j,i] <- max(getFunctionalPPA(prf, logBF))
    }
    if (any(is.na(weightedPRFMat[,i]))) {
      minWeightedPRF <- min(weightedPRFMat[,i], na.rm=T)
      minWeightedMaxPRF <- min(weightedMaxPRFMat[,i], na.rm=T)
      minMaxPRF <- min(maxPRFMat[,i], na.rm=T)
      weightedPRFMat[is.na(weightedPRFMat[,i]),i] <- minWeightedPRF
      weightedMaxPRFMat[is.na(weightedMaxPRFMat[,i]),i] <- minWeightedMaxPRF
      maxPRFMat[is.na(maxPRFMat[,i]),i] <- minMaxPRF
    }
  }
  if (myargs$debug) {
    write.table(weightedPRFMat, "weightedPRFMat.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(weightedMaxPRFMat, "weightedMaxPRFMat.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPRFMat, "maxPRFMat.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPPAMat, "maxPPAMat.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(funcWeightedPRFMat, "funcWeightedPRFMat.txt", row.names=T, col.names=T, quote=F, sep="\t")
  }
  
  weightedPRFMat.sds <- apply(weightedPRFMat, 2, sd)
  #weightedPRFMat.subsetlow <- weightedPRFMat[,weightedPRFMat.sds < median(weightedPRFMat.sds)]
  #weightedPRFMat.subsethigh <- weightedPRFMat[,weightedPRFMat.sds >= median(weightedPRFMat.sds)]
  
  # For each epigenome, get the functional PPA
  weightedPRFMat.rankorder <- getMatrixRankOrdering(weightedPRFMat)
  weightedMaxPRFMat.rankorder <- getMatrixRankOrdering(weightedMaxPRFMat)
  maxPRFMat.rankorder <- getMatrixRankOrdering(maxPRFMat)
  maxPPAMat.rankorder <- getMatrixRankOrdering(maxPPAMat)
  funcWeightedPRFMat.rankorder <- getMatrixRankOrdering(funcWeightedPRFMat)
  #weightedPRFMat.subsetlow.rankorder <- getMatrixRankOrdering(weightedPRFMat.subsetlow)
  #weightedPRFMat.subsethigh.rankorder <- getMatrixRankOrdering(weightedPRFMat.subsethigh)
  if (myargs$debug) {
    write.table(weightedPRFMat.rankorder, "weightedPRFMat.rankorder.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(weightedMaxPRFMat.rankorder, "weightedMaxPRFMat.rankorder.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPRFMat.rankorder, "maxPRFMat.rankorder.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPPAMat.rankorder, "maxPPAMat.rankorder.txt", row.names=T, col.names=T, quote=F, sep="\t")
    write.table(funcWeightedPRFMat.rankorder, "funcWeightedPRFMat.rankorder.txt", row.names=T, col.names=T, quote=F, sep="\t")
    #write.table(weightedPRFMat.subsetlow.rankorder, "weightedPRFMat.subsetlow.rankorder.txt", row.names=T, col.names=T, quote=F, sep="\t")
    #write.table(weightedPRFMat.subsethigh.rankorder, "weightedPRFMat.subsethigh.rankorder.txt", row.names=T, col.names=T, quote=F, sep="\t")
  }
  
  weightedPRF.rra <- aggregateRanks(matrixToCols(weightedPRFMat.rankorder))
  weightedMaxPRF.rra <- aggregateRanks(matrixToCols(weightedMaxPRFMat.rankorder))
  maxPRF.rra <- aggregateRanks(matrixToCols(maxPRFMat.rankorder))
  maxPPA.rra <- aggregateRanks(matrixToCols(maxPPAMat.rankorder))
  funcWeightedPRFMat.rra <- aggregateRanks(matrixToCols(funcWeightedPRFMat.rankorder))
  #weightedPRFMat.subsetlow.rra <- aggregateRanks(matrixToCols(weightedPRFMat.subsetlow.rankorder))
  #weightedPRFMat.subsethigh.rra <- aggregateRanks(matrixToCols(weightedPRFMat.subsethigh.rankorder))
  
  weightedPRFAndMaxPRF.rra <- aggregateRanks(matrixToCols(cbind(weightedPRFMat.rankorder, maxPRFMat.rankorder)))
  
  if (myargs$debug) {
    print("Weighted PRF sum score ranking")
    write.table(weightedPRF.rra, file="", row.names=F, quote=F, sep="\t")
    print("Weighted max PRF ranking")
    write.table(weightedMaxPRF.rra, file="", row.names=F, quote=F, sep="\t")
    print("Max PRF score ranking")
    write.table(maxPRF.rra, file="", row.names=F, quote=F, sep="\t")
    print("Max functional PPA ranking")
    write.table(maxPPA.rra, file="", row.names=F, quote=F, sep="\t")
    print("FuncPPA weighted PRF sum score ranking")
    write.table(funcWeightedPRFMat.rra, file="", row.names=F, quote=F, sep="\t")
  }
  
  colnames(weightedPRF.rra) <- c("Name", "Weighted PRF sum score")
  colnames(weightedMaxPRF.rra) <- c("Name", "Weighted max PRF")
  colnames(maxPRF.rra) <- c("Name", "Max PRF score")
  colnames(maxPPA.rra) <- c("Name", "Max functional PPA")
  colnames(funcWeightedPRFMat.rra) <- c("Name", "FuncPPA weighted PRF sum score")
  
  weightedPRF.eidNames <- metadata.df[weightedPRF.rra$Name, "Epigenome.Mnemonic"]
  weightedMaxPRF.eidNames <- metadata.df[weightedMaxPRF.rra$Name, "Epigenome.Mnemonic"]
  maxPRF.eidNames <- metadata.df[maxPRF.rra$Name, "Epigenome.Mnemonic"]
  maxPPA.eidNames <- metadata.df[maxPPA.rra$Name, "Epigenome.Mnemonic"]
  funcWeightedPRFMat.eidNames <- metadata.df[funcWeightedPRFMat.rra$Name, "Epigenome.Mnemonic"]
  
  weightedPRF.rra$Name <- paste0(weightedPRF.rra$Name, ", ", weightedPRF.eidNames)
  weightedMaxPRF.rra$Name <- paste0(weightedMaxPRF.rra$Name, ", ", weightedMaxPRF.eidNames)
  maxPRF.rra$Name <- paste0(maxPRF.rra$Name, ", ", maxPRF.eidNames)
  maxPPA.rra$Name <- paste0(maxPPA.rra$Name, ", ", maxPPA.eidNames)
  funcWeightedPRFMat.rra$Name <- paste0(funcWeightedPRFMat.rra$Name, ", ", funcWeightedPRFMat.eidNames)
  
  combined.rra <- cbind(weightedPRF.rra, weightedMaxPRF.rra, maxPRF.rra, maxPPA.rra, funcWeightedPRFMat.rra)
  print("Combined rankings table")
  write.table(combined.rra, file="", row.names=F, quote=F, sep="\t")

  if (myargs$perlocus) {
    locus.df <- as.data.frame(maxPRFMat)
    for (locus in colnames(locus.df)) {
      prfOrderDecr <- order(locus.df[,locus], decreasing=T)
      ordered_prfs <- locus.df[prfOrderDecr, locus]
      ordered_eidNames <- metadata.df[rownames(locus.df)[prfOrderDecr], "Epigenome.Mnemonic"]
      ordered_eids <- metadata.df[rownames(locus.df)[prfOrderDecr], "Epigenome.ID"]
      locus.df[,locus] <- paste0(ordered_prfs, ", ", ordered_eidNames, ", ", ordered_eids)
    }
    print("")
    print("Locus rankings table")
    write.table(locus.df, file="", row.names=F, quote=F, sep="\t")
  }
}


getMatrixRankOrdering = function(mat)
{
  nrows <- dim(mat)[1]
  ncols <- dim(mat)[2]
  rankOrderMatrix <- matrix(nrow=nrows, ncol=ncols)
  for (i in 1:ncols) {
    rankOrderMatrix[,i] <- names(mat[order(mat[,i],decreasing=T), i])
  }
  colnames(rankOrderMatrix) <- colnames(mat)
  rankOrderMatrix
}

matrixToCols = function(mat)
{
  lapply(seq_len(ncol(mat)), function(i) mat[,i])
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

getBFFromPval = function(df)
{
  df[,'Z'] <- pToZStat(df$Pval)
  df$F[is.na(df$F)] <- 0.01
  df$N[is.na(df$N)] <- 1000
  apply(df[,c('Z','F','N'),with = FALSE], 1, function(d) calcLogBF(d['Z'],d['F'],d['N']))
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

getSegmentNaivePPAs = function(df)
{
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

sumlog = function(logx, logy)
{
  if (logx > logy) return(logx + log(1 + exp(logy-logx)))
  else return(logy + log(1 + exp(logx-logy)))
}


getSegmentFunctionalPPAs = function(df, PRFcol)
{
  segIndices <- getUniqueIndices(df$segment)
  ppas <- vector(mode="double")
  # For each gene/segment...
  for (i in 1:length(segIndices)) {
    startIndex <- segIndices[[i]][[2]]
    endIndex <- segIndices[[i]][[3]]
    ppas[startIndex:endIndex] <- getFunctionalPPA(df[startIndex:endIndex,PRFcol,with=F], df$logBF[startIndex:endIndex])
  }
  ppas
}

getFunctionalPPA = function(vecPRF, vecLogBF)
{
  sumLogPRF <- -1000
  for (prfScore in vecPRF) {
    if (!is.na(prfScore)) {
      sumLogPRF <- sumlog(sumLogPRF, prfScore)
    }
  }
  vecSnpPri <- vecPRF - sumLogPRF
  vecSnpPriBFSum <- vecSnpPri + vecLogBF
  
  logsegbf <- -1000
  for (val in vecSnpPriBFSum) {
    if (!is.na(val)) {
      logsegbf <- sumlog(logsegbf, val)
    }
  }
  
  # Assume fine mapping, so p(assoc) in segment is 1
  segPPA <- 1
  vecPPA <- exp(vecSnpPriBFSum - logsegbf)
  vecPPA[is.na(vecPPA)] <- 0
  vecPPA
}



###########################################################################
##' commandArgs parsing
##' 
##' return a named list of command line arguments
##'
##' Usage:
##' call the R script thus
##'   ./myfile.R --args myarg=something
##' or
##'   R CMD BATCH --args myarg=something myfile.R
##'
##' Then in R do
##'   myargs <- getArgs()
##' and myargs will be a named list
##' > str(myargs)
##' List of 2
##' $ file : chr "myfile.R"
##' $ myarg: chr "something"
##'
##' @title getArgs
##' @param verbose print verbage to screen 
##' @param defaults a named list of defaults, optional
##' @return a named list
##' @author Chris Wallace
getArgs = function(verbose=FALSE, defaults=NULL) {
  myargs <- gsub("^--","",commandArgs(TRUE))
  setopts <- !grepl("=",myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
  myargs.list <- strsplit(myargs,"=")
  myargs <- lapply(myargs.list,"[[",2 )
  names(myargs) <- lapply(myargs.list, "[[", 1)
  
  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE
  
  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }
  
  ## verbage
  if(verbose) {
    cat("read",length(myargs),"named args:\n")
    print(myargs)
  }
  myargs
}


###########################################################################

main()

