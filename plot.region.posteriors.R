#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(data.table)
library(grid)
library(gridExtra)
library(plyr)

options(datatable.verbose = F)
options(datatable.showProgress = F)
options(datatable.fread.datatable=T)
options(stringsAsFactors = F)

myargs <- NULL
metadata.df <- NULL
numTopSNPsToLabel <- 6
plotRegion <- F
plotZoomin <- F
plotCandidates <- F
prfBarPlot <- T

locusTitleStr <- character(0)
addChrPosToTitle <- T
fineMappingTitleStr <- ""

main = function()
{
  myargs <<- getArgs()
  print(myargs)
  
  if (is.null(myargs$outstem)) {
    stop("Missing required column: 'outstem' argument should give the base path for output files.")
  }
  if (is.null(myargs$file)) {
    #stop("Missing parameter --file")
    setwd("/Users/js29/Downloads")
    myargs$file <<- "IBD_finemapping.PRF.txt"
  }
  if (is.null(myargs$verbose)) {
    myargs$verbose <<- F
  }
  
  df <- fread(myargs$file, header=T, sep='\t', na.strings=c("NA","N/A",""))
  
  if (is.null(myargs$segmentCol)) {
    stop("Missing required column: 'segmentCol' argument should give the name of the region/segment column.")
  }
  if (is.null(myargs$maxSegments)) {
    myargs$maxSegments <<- 1e6
  }
  if (is.null(myargs$pCol)) {
    stop("Missing required column: 'pCol' argument should give the name of the P value column.")
  }
  if (is.null(myargs$F)) {
    stop("Missing required column: 'F' argument should give the name of the allele frequency column.")
  }
  if (!is.null(myargs$NCol)) {
    df$N = as.integer(df[,myargs$NCol,with=F][[1]])
  } else if (!is.null(myargs$N)) {
    df$N = as.integer(myargs$N)
  } else if (! ("N" %in% colnames(df))) {
    stop("Missing required column: either column 'N' should be present, or 'N' argument should give the sample size N, or 'NCol' argument should give the name of the sample size column.")
  }
  if (!is.null(myargs$debug)) {
    myargs$debug <<- T
  } else {
    myargs$debug <<- F
  }
  if (is.null(myargs$scoreCol)) {
    stop("Missing required column: 'scoreCol' argument should give the name of the PRF score column.")
  }
  if (is.null(myargs$snpCol)) {
    stop("Missing required column: 'snpCol' argument should give the name of the snp rsID column.")
  }
  if (is.null(myargs$chrCol)) {
    stop("Missing required column: 'chrCol' argument should give the name of the chromosome column.")
  }
  if (is.null(myargs$posCol)) {
    stop("Missing required column: 'posCol' argument should give the name of the chromosomal position column.")
  }
  if (is.null(myargs$metadata)) {
    myargs$metadata = "/lustre/scratch110/sanger/dg13/share/jeremy/prfscores/epigenomes/epigenome.metadata.txt"
  }
  if (!is.null(myargs$locusTitle)) {
    locusTitleStr <<- myargs$locusTitle
    addChrPosToTitle <<- F
  }
  showAllSnps = F
  if (!is.null(myargs$showAllSnps)) {
    showAllSnps = T
  }
  print(myargs)
  metadata.df <<- read.delim(myargs$metadata, stringsAsFactors=F, header=T)
  metadata.df <<- metadata.df[order(metadata.df$Epigenome.ID),]
  rownames(metadata.df) <<- metadata.df$Epigenome.ID
  
  df$segment <- df[,myargs$segmentCol,with=F][[1]]
  df$Pval <- as.double(df[,myargs$pCol,with=F][[1]])
  df$Pval[df$Pval == 0] <- 1e-300
  df$F <- as.double(df[,myargs$F,with=F][[1]])
  df$score <- as.double(df[,myargs$scoreCol,with=F][[1]])
  df$snpid <- df[,myargs$snpCol,with=F][[1]]
  df$chr <- df[,myargs$chrCol,with=F][[1]]
  df$pos <- as.integer(df[,myargs$posCol,with=F][[1]])
  
  # Get approximate Bayes Factors for each SNP
  df[,"logBF"] <- getBFFromPval(df)
  
  # Get the posterior probability of association for each SNP
  # Requires logBF to have been defined
  df[,"PPANaive"] <- getSegmentNaivePPAs(df)
  if (!showAllSnps) {
    print("Subsetting to SNPs with Bayes factor >= 0")
    df <- df[df$logBF >= 0,]
  }
  
  segNames <- unique(df$segment)
  segIndices <- getUniqueIndices(df$segment)
  
  pdf(file=paste0(myargs$outstem, ".pdf"), width=12, height=10)
  numSegs <- min(length(segIndices), as.integer(myargs$maxSegments))
  for (i in 1:numSegs) {
    startIndex <- segIndices[[i]][[2]]
    endIndex <- segIndices[[i]][[3]]
    numSnps <- endIndex - startIndex + 1
    segname <- df[startIndex,]$segment
    df.seg <- df[startIndex:endIndex,]
    df.seg <- df.seg[order(df.seg$logBF, decreasing=T)]
    df.seg$gwasRank <- as.character(1:numSnps)
    
    #df.seg$pi <- getFunctionalPrior(df.seg$score)
    df.seg$PPA <- getFunctionalPPA(df.seg$score, df.seg$logBF)
    
    print(paste("Segment:", segname))
    eid <- gsub("_PRF", "", myargs$scoreCol)
    eidName <- metadata.df[eid, "Roadmap.Epigenome.name"]
    fineMappingTitleStr <<- paste0("\nFine mapping with ", eid, "-", eidName)
      
    if (plotRegion) {
      if (length(locusTitleStr) <= 0) {
        #locusTitleStr <<- paste0("FULL LOCUS: ", segname)
        locusTitleStr <<- "FULL LOCUS: "
      }
      makePlot(df.seg, TRUE)
    }
    
    # Get SNPs within a logBF of 5 of the lead (but must be >= 2)
    # First find the index of the lead SNP
    leadIdx <- which(df.seg$logBF == max(df.seg$logBF))
    leadLogBF <- df.seg[leadIdx]$logBF
    
    # Get SNPs in this segment, and consider only those in the 95% statistical credible set
    df.seg.cred <- df.seg[order(PPANaive, decreasing=T),]
    if (!showAllSnps) {
      df.seg.cred <- df.seg.cred[PPANaive > 0.001,]
    }
    df.seg.cred <- df.seg.cred[, cumPPANaive := cumsum(PPANaive)]
    print(dim(df.seg.cred))
    credSetInd <- which(df.seg.cred$cumPPANaive <= 0.95)
    if (length(credSetInd) == 0) {
      maxCredSetInd <- 1
    } else {
      # Be sure to include the SNP that puts us just past the 95% credible set level
      maxCredSetInd = max(credSetInd)
      if (max(credSetInd) < dim(df.seg.cred)[1]) {
        maxCredSetInd = maxCredSetInd + 1
      }
    }
    # For saving the credible set, be strict and include only SNPs actually in the credible set
    # Recalculate the PPAs using only these SNPs
    df.seg.cred.save <- df.seg.cred[1:maxCredSetInd]
    df.seg.cred.save$PPA <- getFunctionalPPA(df.seg.cred.save$score, df.seg.cred.save$logBF)
    df.seg.cred.o <- df.seg.cred.save[order(-PPA)]
    df.seg.cred.save <- df.seg.cred.o[, cumPPA := cumsum(PPA)][order(as.integer(gwasRank))]
    # Write a table with all SNPs in the credible set
    write.table(df.seg.cred.save[1:maxCredSetInd], paste0(myargs$outstem, ".credibleset.95.txt"), quote=F, row.names=F, sep="\t", col.names=T)

    # When plotting include at least the top 5 SNPs
    minToInclude <- min(5, dim(df.seg.cred)[1])
    credSetMaxIndex <- max(maxCredSetInd, minToInclude)

    df.seg.cred <- df.seg.cred[1:credSetMaxIndex]
    df.seg.cred$PPA <- getFunctionalPPA(df.seg.cred$score, df.seg.cred$logBF)
    df.seg.cred.o <- df.seg.cred[order(-PPA)]
    df.seg.cred <- df.seg.cred.o[, cumPPA := cumsum(PPA)][order(as.integer(gwasRank))]
    
    credSetMinPos <- min(df.seg.cred$pos)
    credSetMaxPos <- max(df.seg.cred$pos)
    
    if (plotZoomin) {
      df.zoomin <- df.seg[pos >= credSetMinPos & pos <= credSetMaxPos]
      if (length(locusTitleStr) <= 0) {
        #locusTitleStr <<- paste0("CREDIBLE SET INTERVAL: ", segname)
        locusTitleStr <<- "CREDIBLE SET INTERVAL: "
      }
      makePlot(df.zoomin, TRUE)
    }
    if (plotCandidates) {
      if (length(locusTitleStr) <= 0) {
        #locusTitleStr <<- paste0("CREDIBLE SET SNPS ONLY: ", segname)
        locusTitleStr <<- "CREDIBLE SET SNPS ONLY: "
      }
      makePlot(df.seg.cred, TRUE)
    }
    
    if (prfBarPlot) {
      makePRFBarPlot(df.seg.cred)
    }

    ## Doing the same as above with facets
    #makeFacetPlot(df.candidate)  
  }
  dev.off()
  print(warnings())
}


makePlot = function(plot.df, labelSNPText)
{
  if (myargs$debug) {
    print("plot.df dimensions:")
    print(dim(plot.df))
    print(paste("labelSNPText:", labelSNPText))
  }
  plot.df$pi <- getFunctionalPrior(plot.df$score)
  plot.df$PPA <- getFunctionalPPA(plot.df$score, plot.df$logBF)

  # Identify SNPs to label -- ones that are high scoring in either the GWAS
  # or the posteriors
  numSnps <- dim(plot.df)[1]
  end <- min(numSnps, numTopSNPsToLabel)
  topGWASSNPs <- order(plot.df$logBF, decreasing=T)[1:end]
  topPPASNPs <- order(plot.df$PPA, decreasing=T)[1:end]
  labelSNPIdxs <- unique(c(topGWASSNPs, topPPASNPs))
  plot.df$ptSize <- 0.5
  plot.df$labelID <- ""
  plot.df[labelSNPIdxs,]$labelID <- plot.df[labelSNPIdxs,]$gwasRank
  plot.df[labelSNPIdxs,]$ptSize <- 1.5
  
  labelSNPs <- plot.df[labelSNPIdxs, ]

  if (addChrPosToTitle) {
    plotTitle <- paste0(locusTitleStr, plot.df$chr[1], ":", min(plot.df$pos), "-", max(plot.df$pos), fineMappingTitleStr)
  } else {
    plotTitle <- paste0(locusTitleStr, fineMappingTitleStr)
  }
  
  angl = 0
  vj = 0.5
  hj = -0.5
  sz = 15
  gwas.plot <- ggplot(plot.df, aes(x=pos, y=logBF, colour=labelID, size=ptSize)) + geom_point() + scale_size(range=c(2,4)) + theme_bw(sz) + theme(plot.title=element_text(size=sz, hjust=0.98, vjust=0.5), axis.text.x=element_blank(), axis.title.x=element_blank(), plot.margin=unit(c(0,3,0,3),"mm"))
  gwas.plot <- gwas.plot + geom_point(data=labelSNPs, aes(colour=labelID, size=4))
  gwas.plot <- gwas.plot + ggtitle("GWAS") + ylab("ln(BF)") + theme(legend.position="none")
  #gwas.plot <- ggplot(plot.df, aes(x=pos, y=logBF, colour=labelID, size=ptSize)) + geom_point() + scale_size(range=c(2,3)) + theme_grey(sz) + theme(plot.title=element_text(size=sz, hjust=0.98, vjust=-2.5), axis.text.x=element_blank(), axis.title.x=element_blank(), plot.margin=unit(c(0,3,0,3),"mm"))
  #gwas.plot <- gwas.plot + ggtitle("RA GWAS") + ylab("ln(BF)") + geom_dl(aes(label=labelID),list("smart.grid", cex = 0.8)) + theme(legend.position="none") 
  
  gwas.naive.plot <- ggplot(plot.df, aes(x=pos, y=PPANaive, colour=labelID, size=ptSize)) + geom_point() + scale_size(range=c(2,4)) + theme_bw(sz) + theme(plot.title=element_text(size=sz, hjust=0.98, vjust=0.5), axis.text.x=element_blank(), axis.title.x=element_blank(), plot.margin=unit(c(0,3,0,3),"mm"))
  gwas.naive.plot <- gwas.naive.plot + geom_point(data=labelSNPs, aes(colour=labelID, size=4))
  gwas.naive.plot <- gwas.naive.plot + ggtitle("Naive GWAS Posterior") + ylab("Naive PPA") + theme(legend.position="none")

  prf.plot <- ggplot(plot.df, aes(x=pos, y=score, colour=labelID, size=ptSize)) + geom_point() + scale_size(range=c(2,4)) + theme_bw(sz) + theme(plot.title = element_text(size=sz, hjust=0.98, vjust=0.5), axis.text.x=element_blank(), axis.title.x=element_blank(), plot.margin=unit(c(0,3,0,3),"mm"))
  prf.plot <- prf.plot + geom_point(data=labelSNPs, aes(colour=labelID, size=4))
  eid <- gsub("_PRF", "", myargs$scoreCol)
  eidName <- metadata.df[eid, "Roadmap.Epigenome.name"]
  prf.plot <- prf.plot + ggtitle(eidName) + ylab("PRF Score") +
              theme(legend.position="none") +
              theme(plot.title = element_text(colour="blue"))
  prfMaxMinLabel <- paste(max(plot.df$score), "/", plot.df$score[labelSNPIdxs[1]])
  prf.plot <- prf.plot + annotate("text", x = max(plot.df$pos), y = max(plot.df$score), label = prfMaxMinLabel, angle=angl, vjust=1, hjust=1, size=4)
  
  ppa.plot <- ggplot(plot.df, aes(x=pos, y=PPA, colour=labelID, size=ptSize)) + geom_point() + scale_size(range=c(2,4)) + theme_bw(sz) + theme(plot.title = element_text(size=sz, hjust=0.98, vjust=0.5), plot.margin=unit(c(0,3,0,3),"mm"))
  ppa.plot <- ppa.plot + geom_point(data=labelSNPs, aes(colour=labelID, size=4))
  ppa.plot <- ppa.plot + ggtitle(paste0("Posterior - ", eidName)) + ylab("PPA") +
              theme(legend.position="none") +
              theme(plot.title = element_text(colour="blue"))

  if (labelSNPText) { #!is.null(labelSNPs)) {
    if (length(labelSNPs$pos) > 0) {
      gwas.plot <- gwas.plot + annotate("text", x = labelSNPs$pos, y = labelSNPs$logBF, label = labelSNPs$labelID, angle=angl, vjust=vj, hjust=hj, size=4)
      prf.plot <- prf.plot + annotate("text", x = labelSNPs$pos, y = labelSNPs$score, label = labelSNPs$labelID, angle=angl, vjust=vj, hjust=hj, size=4)
      ppa.plot <- ppa.plot + annotate("text", x = labelSNPs$pos, y = labelSNPs$PPA, label = labelSNPs$labelID, angle=angl, vjust=vj, hjust=hj, size=4)
      gwas.naive.plot <- gwas.naive.plot + annotate("text", x = labelSNPs$pos, y = labelSNPs$PPANaive, label = labelSNPs$labelID, angle=angl, vjust=vj, hjust=hj, size=4)
    }
  }

  plots <- list(gwas.plot, gwas.naive.plot, prf.plot)
  #plots <- list(prf.plot)
  extraPrfPlots <- list()
  if (!is.null(myargs$extraScores)) {
    scoreNames <- unlist(strsplit(myargs$extraScores, ",", fixed=T))
    for (scoreName in scoreNames) {
      plot.df$score2 <- as.double(unlist(plot.df[,scoreName,with=F]))
      plot.df$pi2 <- getFunctionalPrior(plot.df$score2)
      labelSNPs <- plot.df[labelSNPIdxs,]
      
      prf.plot2 <- ggplot(plot.df, aes(x=pos, y=score2, colour=labelID, size=ptSize)) + geom_point() + scale_size(range=c(2,4)) + theme_bw(sz) + theme(plot.title = element_text(size=sz, hjust=0.98, vjust=0.5), axis.text.x=element_blank(), axis.title.x=element_blank(), plot.margin=unit(c(0,3,0,3),"mm"))
      prf.plot2 <- prf.plot2 + geom_point(data=labelSNPs, aes(colour=labelID, size=4))
      eid <- gsub("_PRF", "", scoreName)
      eidName <- metadata.df[eid, "Roadmap.Epigenome.name"]
      prf.plot2 <- prf.plot2 + ggtitle(eidName) + ylab("PRF score") +
                  theme(legend.position="none")
      prfMaxMinLabel <- paste(max(plot.df$score2), "/", plot.df$score2[labelSNPIdxs[1]])
      prf.plot2 <- prf.plot2 + annotate("text", x = max(plot.df$pos), y = max(plot.df$score2), label = prfMaxMinLabel, angle=angl, vjust=1, hjust=1, size=4)
                                      
      if (labelSNPText) { #!is.null(labelSNPs)) {
        if (length(labelSNPs$pos) > 0) {
          prf.plot2 <- prf.plot2 + annotate("text", x = labelSNPs$pos, y = labelSNPs$score2, label = labelSNPs$labelID, angle=angl, vjust=vj, hjust=hj, size=4)
        }
      }
      
      extraPrfPlots <- append(extraPrfPlots, list(prf.plot2))
      plots <- append(plots, list(prf.plot2))
    }
  }
  plots <- append(plots, list(ppa.plot))
  grid.arrange(grobs=plots, ncol=1, top=textGrob(plotTitle,gp=gpar(fontface="bold")))
}

makePRFBarPlot = function(plot.df)
{
  # Identify SNPs to label -- ones that are high scoring in either the GWAS
  # or the posteriors
  numSnps <- dim(plot.df)[1]
  end <- min(numSnps, numTopSNPsToLabel)
  topGWASSNPs <- order(plot.df$logBF, decreasing=T)[1:end]
  topPPASNPs <- order(plot.df$PPA, decreasing=T)[1:end]
  labelSNPIdxs <- unique(c(topGWASSNPs, topPPASNPs))
  plot.df$labelID <- ""
  plot.df[labelSNPIdxs,]$labelID <- plot.df[labelSNPIdxs,]$gwasRank
  bar.df <- plot.df[labelSNPIdxs, ]
  
  if (addChrPosToTitle) {
    plotTitle <- paste0(locusTitleStr, plot.df$chr[1], ":", min(plot.df$pos), "-", max(plot.df$pos), fineMappingTitleStr)
  } else {
    plotTitle <- paste0(locusTitleStr, fineMappingTitleStr)
  }
  
  col.names <- colnames(bar.df)  
  bar.df$snpid <- sapply(bar.df$snpid, FUN=function(x) strsplit(x, ";", fixed=T)[[1]][1])
  
  rank_ids <- paste0(bar.df$gwasRank, ":\n", bar.df$snpid)
  bar.df$rank_id <- factor(rank_ids, levels=rank_ids[order(as.integer(bar.df$gwasRank))])
  if (length(grep("x.", col.names, fixed=T)) <= 0) {
    return()
  }
  prfColnames <- col.names[grep("x.", col.names, fixed=T)]
  bar.df <- bar.df[,c("rank_id",prfColnames),with=F]
  #colnames(bar.df) <- c("snpid","TSS Dist","DNase","UTR3","UTR5","FantomEnh","EncodeEnh","EncodeRepr","H3k4me1","H3k4me3","H3k27ac","H3k27me3","H3k36me3","Gerp","Effect-SNP")
  
  library(reshape2)
  melt.df <- melt(bar.df, id.vars=c("rank_id"), measure.vars=prfColnames)
  melt.df$variable <- gsub("^x.", "", melt.df$variable)
  melt.df$variable <- gsub("DNAMethylSBS.fm", "DNAMethyl", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("intron.samegene", "intron.same", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("intron.diffgene", "intron.diff", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("coding.samegene", "coding.same", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("UTR5.samegene", "UTR5.same", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("UTR3.samegene", "UTR3.same", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("UTR3.diffgene", "UTR3.diff", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("UTR5.diffgene", "UTR5.diff", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("effect-snp.nmotifs", "effect-snp", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("switch-snp.nmotifs", "switch-snp", melt.df$variable, ignore.case=T)
  melt.df$variable <- gsub("Enh.Fantom.tpm", "Enh.Fantom", melt.df$variable, ignore.case=T)
  dat.pos <- subset(melt.df, value > 0)
  dat.neg <- subset(melt.df, value < 0)
  #dat.neg <- dat.neg[variable != "H3k4me3"]
  
  #bar.colors <- read.table("prf.bar.colors.sorted.txt",header=T,sep="\t",comment.char="%",stringsAsFactors=F)
  #bar.colors$variable = factor(bar.colors$variable, levels=c("TSS Dist","DNase","UTR3","UTR5","FantomEnh","EncodeEnh","EncodeRepr","H3k4me1","H3k4me3","H3k27ac","H3k27me3","H3k36me3","Gerp","Effect-SNP"))
  #colorValues <- bar.colors$color
  #names(colorValues) <- bar.colors$variable
  
  g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  dat.pos <- ddply(dat.pos, .(rank_id), transform, pos = cumsum(value) - (0.5 * value))
  dat.neg <- ddply(dat.neg, .(rank_id), transform, pos = cumsum(value) - (0.5 * value))
  #geom_text(aes(label = variable, y = pos), size = 3) +
  
  dat.pos.lab <- dat.pos[dat.pos$value >= 0.2,]
  dat.neg.lab <- dat.neg[dat.neg$value <= -0.2,]
  
  firstVar = dat.pos$rank_id[1]
  stackedPlot <- ggplot() +
    geom_bar(data = dat.pos, aes(x=rank_id, y=value, fill=variable), stat = "identity") +
    geom_bar(data = dat.neg, aes(x=rank_id, y=value, fill=variable), stat = "identity") +
    theme_bw(14) +
    geom_hline(aes(yintercept=0)) +
    theme(legend.position="right", legend.key.height=unit(0.6, "cm")) + 
    ylab("Annotation contributions") + 
    annotate("text", x = dat.pos.lab$rank_id, y=dat.pos.lab$pos, label=dat.pos.lab$variable, size=4) +
    annotate("text", x = dat.neg.lab$rank_id, y=dat.neg.lab$pos, label=dat.neg.lab$variable, size=4)
  #annotate("text", x = dat.pos$rank_id[1], y = dat.pos[dat.pos$rank_id==firstVar&dat.pos$variable=="TSSDist",]$value, label = c("my label"))
    
  stackedPlot <- stackedPlot + theme_bw(15) +
    theme(legend.position="none") +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(plot.margin=unit(c(0,0,0.2,0.1),"cm"))
    
  
  bar.df$total <- apply(as.data.frame(bar.df)[,prfColnames], MARGIN=1, FUN = function(x) {sum(x)})
  sumPlot <- ggplot() + 
    geom_bar(data = bar.df, aes(x=rank_id, y=total), fill="#3366FF", stat = "identity") +
    theme_bw(14) + theme(legend.position="none", axis.title.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=unit(c(0.3,0,0,0.3),"cm")) + ylab("PRF score")
  
  #annotlegend <- g_legend(stackedPlot)
  #grid.arrange(arrangeGrob(sumPlot, stackedPlot, ncol=1, heights=c(2,5), padding=unit(0, "cm")), annotlegend, widths=c(5,1), ncol=2, top=textGrob(plotTitle,gp=gpar(fontface="bold")))
  
  grid.arrange(arrangeGrob(sumPlot, stackedPlot, ncol=1, heights=c(2,5), padding=unit(0, "cm")), ncol=1, top=textGrob(plotTitle,gp=gpar(fontface="bold")))
}

makeFacetPlot = function(plot.df, labelSNPs)
{
  plot.df <- prf.candidate
  facet1 <- data.frame(pos=plot.df$pos, val=plot.df$logBF)
  facet1$Type = "GWAS"
  facet2 <- data.frame(pos=plot.df$pos, val=plot.df$pi)
  facet2$Type = "PRF score"
  facet3 <- data.frame(pos=plot.df$pos, val=plot.df$PPA)
  facet3$Type = "PPA"
  
  prf.facet <- rbind(as.data.frame(facet1), as.data.frame(facet2), as.data.frame(facet3))
  prf.facet <- rbind(facet1, facet2, facet3)
  prf.facet$Type <- factor(prf.facet$Type, levels=c("GWAS", "PRF score", "PPA"))
  facet.plot <- ggplot(prf.facet, aes(x=pos, y=val)) + geom_point() + theme_grey(16) + theme(plot.title=element_text(size=16), axis.text=element_text(colour="black"))
  #facet.plot <- ggplot(prf.facet, aes(x=pos, y=val)) + geom_point()
  facet.plot <- facet.plot + facet_wrap(~Type, ncol=1, scales="free_y")
  print(facet.plot)
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

getFunctionalPrior = function(vecPRF)
{
  sumLogPRF <- -1000
  for (prfScore in vecPRF) {
    if (!is.na(prfScore)) {
      sumLogPRF <- sumlog(sumLogPRF, prfScore)
    }
  }
  vecSnpPri <- vecPRF - sumLogPRF
  vecSnpPri
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

                                                                                                                                                                                                                                                                