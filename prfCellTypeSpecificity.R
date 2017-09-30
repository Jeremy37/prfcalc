#!/usr/bin/env Rscript
library(data.table)
library(RobustRankAggreg)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
#source("ppaFunctions.R")
options(stringsAsFactors = F)

myargs = NULL
epigenomeIDs = NULL
PRFCALC_DIR = "/nfs/teams/team170/jeremy/src/prfcalc/"
PRF_REGRESSION_WINDOW = 10000

main = function()
{
  myargs <<- getArgs()
  print(myargs)
  
  if (is.null(myargs$file)) {
    stop("Missing parameter --file")
  }
  if (is.null(myargs$output)) {
    stop("Missing parameter --output")
  }
  myargs$verbose = asBooleanArg(myargs$verbose)
  myargs$debug = asBooleanArg(myargs$debug)
  
  if (is.null(myargs$metadata)) {
    myargs$metadata = paste0(PRFCALC_DIR, "epigenome.metadata.txt")
  }
  metadata.df = read.delim(myargs$metadata, stringsAsFactors=F, header=T)
  metadata.df = metadata.df[order(metadata.df$Epigenome.ID),]
  metadata.df$eid = metadata.df$Epigenome.ID
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
  myargs$perlocus = asBooleanArg(myargs$perlocus)
  myargs$prfmean = asBooleanArg(myargs$prfmean)
  myargs$prfregression = asBooleanArg(myargs$prfregression)
  myargs$meanadjust = !asBooleanArg(myargs$nomeanadjust)
  myargs$makeplots = !asBooleanArg(myargs$noplots)
  
  if (myargs$prfregression) {
    if (is.null(myargs$posCol)) {
      stop("Missing argument: when --prfregression is used, --posCol must also be specified to indicate the chromosomal position column")
    }
    df$pos = as.integer(df[[myargs$posCol]])
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
    write.table(df, paste0(myargs$output, ".output.temp.txt"), row.names=F, quote=F, sep="\t")
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
  
  if (myargs$meanadjust) {
    eid.prf.df = read.delim(paste0(PRFCALC_DIR, "epigenome.mean.prf.txt"))
    #eid.prf.df = read.delim(paste0("/Users/js29/work/prf2/scoredist/", "epigenome.mean.prf.txt"))
    eid.prf.meanAdjust = eid.prf.df$meanPRF - mean(eid.prf.df$meanPRF)
    names(eid.prf.meanAdjust) = eid.prf.df$eid
  }
  if (myargs$prfregression) {
    df.prfRegression = data.frame()
  }
  
  locusSignificance = double(length = length(segIndices))
  for (i in 1:length(segIndices)) {
    startIndex <- segIndices[[i]][[2]]
    endIndex <- segIndices[[i]][[3]]
    # Get SNPs in this segment, and consider only those in the 95% statistical credible set
    df.seg <- df[startIndex:endIndex,]
    df.seg.cred <- df.seg[PPANaive > 0.01,][order(PPANaive, decreasing=T),]
    df.seg.cred <- df.seg.cred[, cumPPANaive := cumsum(PPANaive)]
    credSetInd <- which(df.seg.cred$cumPPANaive > 0.95)
    if (length(credSetInd) == 0) {
      credSetMax <- nrow(df.seg.cred)
    } else {
      credSetMax <- min(credSetInd)
    }
    df.seg.cred <- df.seg.cred[1:credSetMax]
    if (myargs$verbose) {
      print(sprintf("%s: num SNPs=%d  cred set SNPs=%d", df.seg$segment[1], nrow(df.seg), nrow(df.seg.cred)))
    }
    
    locusSignificance[i] = max(df.seg.cred$logBF)
    
    if (myargs$prfregression) {
      # At each locus we want to take credible set SNPs and other SNPs within
      # a certain distance of the lead SNP, so that they are region-matched.
      leadSNPPos = df.seg.cred$pos[1]
      window.indices = which(abs(df.seg$pos - leadSNPPos) <= PRF_REGRESSION_WINDOW)
      cred.indices = which(df.seg$PPANaive > 0.01)
      indices = unique(c(window.indices, cred.indices))
      if ((length(indices) - length(cred.indices)) < 2) {
        write(paste0("At locus ", df.seg$segment[1], ": too few locus SNPs NOT in credible set - skipping this locus for PRF-PPA regression.\n"), stderr())
      } else {
        df.window = df.seg[indices, ][order(PPANaive, decreasing=T), ]
        df.prfRegression = rbind(df.prfRegression, df.window)
      }
    }
    
    # Get a "baseline" PRF score that respresents some normal PRF score
    # across epigenomes for this region. This is used to scale 
    #prfBaseline = median(unlist(df.seg.cred[,epigenomeIDs,with=F]))
    #df.seg.cred[,epigenomeIDs,with=F] <- df.seg.cred[,epigenomeIDs,with=F] - prfBaseline
    for (j in 1:length(epigenomeIDs)) {
      epigenomeID <- epigenomeIDs[j]
      df.eid <- df.seg.cred[,c("PPANaive", epigenomeID, "logBF"),with=F]
      prfScores <- df.eid[,epigenomeID,with=F][[1]]
      if (myargs$meanadjust) {
        prfScores = prfScores - eid.prf.meanAdjust[epigenomeID]
      }
      
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
  if (myargs$verbose) {
    write.table(weightedPRFMat, paste0(myargs$output, ".weightedPRFMat.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(weightedMaxPRFMat, paste0(myargs$output, ".weightedMaxPRFMat.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPRFMat, paste0(myargs$output, ".maxPRFMat.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPPAMat, paste0(myargs$output, ".maxPPAMat.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(funcWeightedPRFMat, paste0(myargs$output, ".funcWeightedPRFMat.txt"), row.names=T, col.names=T, quote=F, sep="\t")
  }
  
  if (myargs$makeplots) {
    pdf(file=paste0(myargs$output, ".cellTypeSpecificity.plots.pdf"), width=7, height=6)
  }
  
  ## PRF average across loci code
  if (myargs$prfmean) {
    meanMaxPRF = apply(maxPRFMat, 1, mean)
    meanWeightedPRF = apply(weightedPRFMat, 1, mean)
    
    # Also compute the meanMaxPRF and so on after weighting according
    # to each locus' significance (logBF).
    sigWeight = locusSignificance * length(locusSignificance) / sum(locusSignificance)
    meanMaxPRF.bySignificance = apply(maxPRFMat, 1, FUN = function(x) mean(x*sigWeight))
    meanWeightedPRF.bySignificance = apply(weightedPRFMat, 1, FUN = function(x) mean(x*sigWeight))
    prfMeans.df = data.frame(eid = epigenomeIDs,
                             maxPRF = meanMaxPRF,
                             weightedPRF = meanWeightedPRF,
                             maxPRF.bySignificance = meanMaxPRF.bySignificance,
                             weightedPRF.bySignificance = meanWeightedPRF.bySignificance)
    prfMeans.df = prfMeans.df[order(-prfMeans.df$maxPRF),]
    #rownames(prfMeans.df) = epigenomeIDs
    write.table(prfMeans.df, paste0(myargs$output, ".cellTypeSpecificity.avg_across_loci.txt"), row.names=F, col.names=T, quote=F, sep="\t")
    
    if (myargs$makeplots) {
      #plot.df = tidyr::gather(prfMeans.df, key="scoreName", value="scoreVal", c("maxPRF", "weightedPRF", "maxPRF.bySignificance", "weightedPRF.bySignificance"))
      rankScores = function(vec, scoreName) {
        vec = vec[order(-vec)]
        data.frame(eid = names(vec),
                   scoreName = scoreName,
                   score = vec,
                   rank = 1:length(vec))
      }
      plotRankedScores = function(df) {
        tit = df$scoreName[1]
        df = df %>% dplyr::left_join(metadata.df[,c("eid", "Epigenome.shortName", "SuperGroup")], by="eid")
        p = ggplot(df, aes(x=rank, y=score)) + geom_point(aes(color=SuperGroup)) +
          theme_bw(11) + ggtitle(tit) + theme(legend.position="none")
        txt = paste0(sapply(1:6, function(i) {paste0(i, ". ", df$Epigenome.shortName[i])}), collapse="\n")
        p + annotate("text", label = txt, x = max(df$rank) * 0.5, y = max(df$score), vjust = 1, hjust = 0, size=2.5)
      }
      grid.arrange(grobs=list(plotRankedScores(rankScores(meanMaxPRF, "maxPRF")),
                              plotRankedScores(rankScores(meanWeightedPRF, "weightedPRF")),
                              plotRankedScores(rankScores(meanMaxPRF.bySignificance, "maxPRF.bySignificance")),
                              plotRankedScores(rankScores(meanWeightedPRF.bySignificance, "weightedPRF.bySignificance"))),
                   ncol=2)
      if (myargs$verbose) {
        write.table(rankScores(meanMaxPRF, "maxPRF"),
                    file=paste0(myargs$output, ".meanMaxPRF.txt"), row.names=F, col.names=T, quote=F, sep="\t")
        write.table(rankScores(meanWeightedPRF, "weightedPRF"),
                    file=paste0(myargs$output, ".meanWeightedPRF.txt"), row.names=F, col.names=T, quote=F, sep="\t")
      }
      
#       plot.df = rbind(rankScores(meanMaxPRF, "maxPRF"),
#                       rankScores(meanWeightedPRF, "weightedPRF"),
#                       rankScores(meanMaxPRF.bySignificance, "maxPRF.bySignificance"),
#                       rankScores(meanWeightedPRF.bySignificance, "weightedPRF.bySignificance"))
#       p = ggplot(plot.df, aes(x=rank, y=score)) + geom_point() +
#         facet_wrap(~scoreName, scales = "free_y") + theme_bw(14)
    }
  }
  
  ## Regression of PRF vs. PPA for SNPs near association peaks
  if (myargs$prfregression) {
    # At each locus we want to take credible set SNPs and other SNPs within
    # a certain distance of the lead SNP, so that they are region-matched.
    df.reg.summary = data.frame()
    for (j in 1:length(epigenomeIDs)) {
      epigenomeID <- epigenomeIDs[j]
      df.r = data.frame(prf = as.vector(df.prfRegression[[epigenomeID]]),
                        PPA = df.prfRegression$PPANaive)
      s = summary(lm(prf ~ PPA, data=df.r))
      df.reg.summary = rbind(df.reg.summary,
                             data.frame(eid = epigenomeID,
                                        est = s$coefficients["PPA","Estimate"],
                                        pval = s$coefficients["PPA","Pr(>|t|)"]))
    }
    df.reg.summary = df.reg.summary[order(df.reg.summary$pval),]
    df.reg.summary = df.reg.summary %>% dplyr::left_join(metadata.df[,c("eid", "Epigenome.Mnemonic")], by="eid")
    write.table(df.reg.summary, paste0(myargs$output, ".cellTypeSpecificity.prf_regression.txt"), row.names=F, col.names=T, quote=F, sep="\t")
    if (myargs$makeplots) {
      df.reg.summary$log10_pval = -log10(df.reg.summary$pval)
      print(cellTypeRankPlot(df.reg.summary, metadata.df, "log10_pval", "-log10(pval)", "PRF-PPA regression"))
      print(cellTypeRankPlot(df.reg.summary, metadata.df, "est", "slope", "PRF-PPA regression"))
    }
  }
  
  ## Robust rank aggregation code

  # For each epigenome, get the functional PPA
  weightedPRFMat.rankorder <- getMatrixRankOrdering(weightedPRFMat)
  weightedMaxPRFMat.rankorder <- getMatrixRankOrdering(weightedMaxPRFMat)
  maxPRFMat.rankorder <- getMatrixRankOrdering(maxPRFMat)
  maxPPAMat.rankorder <- getMatrixRankOrdering(maxPPAMat)
  funcWeightedPRFMat.rankorder <- getMatrixRankOrdering(funcWeightedPRFMat)
  if (myargs$debug) {
    write.table(weightedPRFMat.rankorder, paste0(myargs$output, ".weightedPRFMat.rankorder.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(weightedMaxPRFMat.rankorder, paste0(myargs$output, ".weightedMaxPRFMat.rankorder.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPRFMat.rankorder, paste0(myargs$output, ".maxPRFMat.rankorder.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(maxPPAMat.rankorder, paste0(myargs$output, ".maxPPAMat.rankorder.txt"), row.names=T, col.names=T, quote=F, sep="\t")
    write.table(funcWeightedPRFMat.rankorder, paste0(myargs$output, ".funcWeightedPRFMat.rankorder.txt"), row.names=T, col.names=T, quote=F, sep="\t")
  }
  
  weightedPRF.rra <- aggregateRanks(matrixToCols(weightedPRFMat.rankorder))
  weightedMaxPRF.rra <- aggregateRanks(matrixToCols(weightedMaxPRFMat.rankorder))
  maxPRF.rra <- aggregateRanks(matrixToCols(maxPRFMat.rankorder))
  maxPPA.rra <- aggregateRanks(matrixToCols(maxPPAMat.rankorder))
  funcWeightedPRFMat.rra <- aggregateRanks(matrixToCols(funcWeightedPRFMat.rankorder))

  weightedPRFAndMaxPRF.rra <- aggregateRanks(matrixToCols(cbind(weightedPRFMat.rankorder, maxPRFMat.rankorder)))
  
  if (myargs$debug) {
    print("Weighted PRF sum score ranking")
    write.table(weightedPRF.rra, file="", row.names=F, quote=F, sep="\t")
    print("Weighted max PRF ranking")
    write.table(weightedMaxPRF.rra, file="", row.names=F, quote=F, sep="\t")
    print("Max PRF score ranking")
    write.table(maxPRF.rra, file="", row.names=F, quote=F, sep="\t")
    # Max functional PPA seems a very poor ranking.
    #print("Max functional PPA ranking")
    #write.table(maxPPA.rra, file="", row.names=F, quote=F, sep="\t")
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
  
  #combined.rra <- cbind(weightedPRF.rra, weightedMaxPRF.rra, maxPRF.rra, maxPPA.rra, funcWeightedPRFMat.rra)
  combined.rra <- cbind(weightedPRF.rra, weightedMaxPRF.rra, maxPRF.rra, funcWeightedPRFMat.rra)
  if (myargs$verbose) {
    print("Combined rankings table")
    write.table(combined.rra, file="", row.names=F, quote=F, sep="\t")
  }
  write.table(combined.rra, file=paste0(myargs$output, ".cellTypeSpecificity.rra.txt"), row.names=F, quote=F, sep="\t")

  if (myargs$makeplots) {
    # Plot of RRA p values across epigenomes
    weightedPRF.rra$eid = rownames(weightedPRF.rra)
    weightedPRF.rra$pval = weightedPRF.rra$`Weighted PRF sum score`
    weightedPRF.rra$log10_pval = -log10(weightedPRF.rra$pval)
    p = cellTypeRankPlot(weightedPRF.rra, metadata.df, "log10_pval", "-log10(pval)", "Robust rank aggregation - weighted PRF")
    print(p)
    
    maxPRF.rra$eid = rownames(maxPRF.rra)
    maxPRF.rra$pval = maxPRF.rra$`Max PRF score`
    maxPRF.rra$log10_pval = -log10(maxPRF.rra$pval)
    p = cellTypeRankPlot(maxPRF.rra, metadata.df, "log10_pval", "-log10(pval)", "Robust rank aggregation - max PRF")
    print(p)
  }
  
  if (myargs$perlocus) {
    doOutputPerLocus = function(mat, scoreName) {
      locus.df <- as.data.frame(mat)
      for (locus in colnames(locus.df)) {
        prfOrderDecr <- order(locus.df[,locus], decreasing=T)
        ordered_prfs <- locus.df[prfOrderDecr, locus]
        ordered_eidNames <- metadata.df[rownames(locus.df)[prfOrderDecr], "Epigenome.Mnemonic"]
        ordered_eids <- metadata.df[rownames(locus.df)[prfOrderDecr], "Epigenome.ID"]
        locus.df[,locus] <- paste0(ordered_prfs, ", ", ordered_eidNames, ", ", ordered_eids)
      }
      if (myargs$verbose) {
        print("")
        print(paste0("Locus rankings table: ", scoreName))
        write.table(locus.df, file="", row.names=F, quote=F, sep="\t")
      }
      write.table(locus.df, file=paste0(myargs$output, ".perlocus.", scoreName, ".maxPRF.txt"), row.names=F, quote=F, sep="\t")
      
      if (myargs$makeplots) {
        # Plot of epigenome ranks across loci
        nloci = ncol(weightedPRFMat.rankorder)
        weightedPRF.perlocus.ranks = as.data.frame(weightedPRFMat.rankorder)
        colnames(weightedPRF.perlocus.ranks) = 1:nloci
        weightedPRF.perlocus.ranks$rank = 1:nrow(weightedPRF.perlocus.ranks)
        weightedPRF.perlocus.ranks = weightedPRF.perlocus.ranks %>% tidyr::gather(key=locus, value=eid, 1:nloci)
        weightedPRF.perlocus.ranks$locus = as.integer(weightedPRF.perlocus.ranks$locus)
        ggplot(weightedPRF.perlocus.ranks, aes(x=locus, y=rank, color=eid)) + geom_line(aes(group = eid)) +
          theme(legend.position="none") +
          ggtitle(paste0("Epigenome ranks per locus: ", scoreName))
        
        library(pheatmap)
        pheatmap(scale(maxPRFMat), show_colnames = T, show_rownames = T, fontsize_row = 8)
        x = hclust(dist(t(maxPRFMat)))
        weightedPRF.perlocus.ranks$locusOrder = weightedPRF.perlocus.ranks$locus
        weightedPRF.perlocus.ranks$locusOrder = as.integer(x$order)
        weightedPRF.perlocus.ranks$eidGroup = "other"
        weightedPRF.perlocus.ranks$eidGroup[weightedPRF.perlocus.ranks$eid == "E066"] = "E066"
        weightedPRF.perlocus.ranks$eidGroup[weightedPRF.perlocus.ranks$eid == "E017"] = "E017"
        weightedPRF.perlocus.ranks$eidGroup[weightedPRF.perlocus.ranks$eid == "E030"] = "E030"
        weightedPRF.perlocus.ranks$eidGroup[weightedPRF.perlocus.ranks$eid == "E124"] = "E124"
        weightedPRF.perlocus.ranks$eidWeight = 1
        weightedPRF.perlocus.ranks$eidWeight[weightedPRF.perlocus.ranks$eid == "E066"] = 3
        weightedPRF.perlocus.ranks$eidWeight[weightedPRF.perlocus.ranks$eid == "E017"] = 1.5
        weightedPRF.perlocus.ranks$eidWeight[weightedPRF.perlocus.ranks$eid == "E030"] = 1.5
        weightedPRF.perlocus.ranks$eidWeight[weightedPRF.perlocus.ranks$eid == "E124"] = 1.5
        ggplot(weightedPRF.perlocus.ranks, aes(x=locusOrder, y=rank, color=eidGroup)) + geom_line(aes(group = eid,size=eidWeight), alpha=0.5) +
          ggtitle(paste0("Epigenome ranks per locus: ", scoreName))
        ggplot(weightedPRF.perlocus.ranks, aes(x=locusOrder, y=rank)) +
          geom_line(aes(group = eid,size=eidWeight), alpha=0.5) +
          geom_line(data=weightedPRF.perlocus.ranks[!weightedPRF.perlocus.ranks$eidGroup %in% c("other", "E066"),],
                    aes(group=eid, color=eidGroup, size=eidWeight)) +
          geom_line(data=weightedPRF.perlocus.ranks[weightedPRF.perlocus.ranks$eidGroup == "E066",],
                    aes(group=eid, color=eidGroup, size=eidWeight)) +
          scale_size(range=c(0.5,3)) +
          scale_y_reverse() + theme_bw() +
          guides(size=FALSE) +
          ggtitle(paste0("Epigenome ranks per locus: ", scoreName))
      }
    }
    doOutputPerLocus(maxPRFMat, "maxPRF")
  }
  
  if (myargs$makeplots) {
    dev.off()
  }
}

cellTypeRankPlot = function(df, metadata.df, scoreCol, scoreName, title)
{
  df$sig = df$pval < 0.05
  df$score = df[,scoreCol]
  df = df[order(-df$score),]
  df$rank = 1:nrow(df)
  if (!("Epigenome.shortName" %in% colnames(df)) & ("eid" %in% colnames(df))) {
    df = df %>% dplyr::left_join(metadata.df[,c("eid", "Epigenome.shortName", "SuperGroup")], by="eid")
  }
  txt = paste0(sapply(1:6, function(i) {paste0(i, ". ", df$Epigenome.shortName[i])}), collapse="\n")
  ggplot(df, aes(x=rank, y=score)) + geom_bar(stat="identity", mapping = aes(fill=SuperGroup)) +
    theme_bw(11) + scale_fill_discrete(name="Tissue group") +
    ylab(scoreName) + ggtitle(title) +
    annotate("text", label = txt, x = max(df$rank) * 0.5, y = max(df$score), vjust = 1, hjust = 0, size=2.5)
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

asBooleanArg = function(arg) {
  if (!is.null(arg)) {
    arg = T
  } else {
    arg = F
  }
  arg
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

