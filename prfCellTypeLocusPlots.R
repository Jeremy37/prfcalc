#!/usr/bin/env Rscript
library(RobustRankAggreg)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(gridExtra)
library(grid)

options(stringsAsFactors = F)
myargs = NULL
PRFCALC_DIR = "/nfs/teams/team170/jeremy/src/prfcalc/"

main = function()
{
  myargs <<- getArgs()
  print(myargs)

  if (is.null(myargs$file)) {
    # This should be the .maxPRFMat.txt or .weightedPRFMat.txt files output
    # by prfCellTypeSpecificity.R
    stop("Missing parameter --file")
  }
  if (is.null(myargs$output)) {
    stop("Missing parameter --output")
  }
  if (is.null(myargs$eid)) {
    stop("Missing parameter --eid")
  }
  if (is.null(myargs$width)) {
    myargs$width = 7
  }
  if (is.null(myargs$height)) {
    myargs$height = 6
  }
  if (!is.null(myargs$detailed)) {
    myargs$detailed = T
  } else {
    myargs$detailed = F
  }
  if (is.null(myargs$metadata)) {
    myargs$metadata = paste0(PRFCALC_DIR, "epigenome.metadata.txt")
  }
  myargs$width = as.numeric(myargs$width)
  myargs$height = as.numeric(myargs$height)
  
  df = read.delim(myargs$file)
  
  eids = strsplit(myargs$eid, split=",", fixed=T)[[1]]
  if (any(!eids %in% rownames(df))) {
    stop(paste0("Not all EIDs specified with --eid (", myargs$eid, ") were found in rownames of input file."))
  }
  
  mat.rankorder <- getMatrixRankOrdering(as.matrix(df))

  pdf(file=paste0(myargs$output, ".heatmap.pdf"), width=9, height=9)
  
  pheatmap(scale(df), show_colnames = T, show_rownames = T, fontsize_row = 4, main = "Clustered heatmap: scaled PRF scores")
  #x = hclust(dist(t(df)))
  #x = hclust(dist(df))
  #plot(x)
  dev.off()
  
  # Plot of epigenome ranks across loci
  nloci = ncol(mat.rankorder)
  perlocus.ranks = as.data.frame(mat.rankorder)
  #colnames(perlocus.ranks) = 1:nloci
  perlocus.ranks$rank = 1:nrow(perlocus.ranks)
  perlocus.ranks = perlocus.ranks %>% tidyr::gather(key=locus, value=eid, 1:nloci)
  #perlocus.ranks$locus = as.integer(perlocus.ranks$locus)

  perlocus.ranks$eidGroup = "other"
  perlocus.ranks$eidWeight = 1
  perlocus.ranks$eidGroup[perlocus.ranks$eid == eids[1]] = eids[1]
  perlocus.ranks$eidWeight[perlocus.ranks$eid == eids[1]] = 3
  for (eid in eids[2:length(eids)]) {
    perlocus.ranks$eidGroup[perlocus.ranks$eid == eid] = eid
    perlocus.ranks$eidWeight[perlocus.ranks$eid == eid] = 1.5
  }
  perlocus.ranks$eidGroup = factor(perlocus.ranks$eidGroup, levels=c(eids, "other"))
  
  pdf(file=paste0(myargs$output, ".loci.pdf"), width=myargs$width, height=myargs$height)
  
  ggplot(perlocus.ranks, aes(x=locus, y=rank)) +
    geom_line(aes(group=eid, size=eidWeight), alpha=0.2) +
    geom_line(data=perlocus.ranks[!perlocus.ranks$eidGroup %in% c("other", eids[1]),],
              aes(group=eid, color=eidGroup, size=eidWeight)) +
    geom_line(data=perlocus.ranks[perlocus.ranks$eidGroup == eids[1],],
              aes(group=eid, color=eidGroup, size=eidWeight)) +
    scale_size(range=c(0.5,3)) +
    scale_y_reverse() + theme_bw() +
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_color_discrete(name="EID") +
    ggtitle("Epigenome ranks per locus")
  
  # DO the same plot but with a log scale on the y (rank) axis
  ggplot(perlocus.ranks, aes(x=locus, y=log2(rank))) +
    geom_line(aes(group=eid, size=eidWeight), alpha=0.2) +
    geom_line(data=perlocus.ranks[!perlocus.ranks$eidGroup %in% c("other", eids[1]),],
              aes(group=eid, color=eidGroup, size=eidWeight)) +
    geom_line(data=perlocus.ranks[perlocus.ranks$eidGroup == eids[1],],
              aes(group=eid, color=eidGroup, size=eidWeight)) +
    scale_size(range=c(0.5,3)) +
    scale_y_reverse() + theme_bw() +
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Epigenome ranks per locus")
  
  epigenome.meta = read.delim(myargs$metadata, stringsAsFactors = F)
  plotsPerPage = 10
  numplots = ceiling(ncol(df) / plotsPerPage)
  for (i in 1:numplots) {
    starti = (i-1)*plotsPerPage + 1
    endi = min(i*plotsPerPage, ncol(df))
    subset.df = data.frame(eid=rownames(df), df[,starti:endi])
    plot.df = tidyr::gather(subset.df, key=locus, value=score, 2:ncol(subset.df))
    
    # Add in epigenome names rather than EIDs
    plot.df = plot.df %>% dplyr::left_join(epigenome.meta %>% rename("eid"="Epigenome.ID") %>% dplyr::select(eid, Epigenome.shortName), by="eid")
    plot.df$locus = gsub("^X", "", plot.df$locus)
    highlighted.eid.names = as.character(unique(plot.df$Epigenome.shortName[plot.df$eid %in% eids]))
    plot.df$Epigenome.shortName[!plot.df$eid %in% eids] = "other"
    plot.df$Epigenome.shortName = factor(plot.df$Epigenome.shortName, levels=c(highlighted.eid.names, "other"))
    plot.df$weight = 1
    plot.df$weight[plot.df$eid %in% eids] = 2
    plot.df$layer0alpha = 2.00 - plot.df$weight
    p = ggplot(plot.df, aes(x=locus, y=score, color=Epigenome.shortName, size=weight)) +
      geom_jitter(aes(size=1, alpha=layer0alpha), width=0.2) +
      geom_jitter(aes(x=locus, y=score, color=Epigenome.shortName, size=weight),
                  data=plot.df %>% dplyr::filter(Epigenome.shortName != "other"),
                  width=0.2, alpha=1) +
      scale_size(range=c(1,4)) + guides(size=FALSE, alpha=FALSE) +
      theme_bw(11) + scale_color_discrete(name="EID") +
      theme(axis.text.x = element_text(angle=45, hjust=1), legend.title = element_blank()) +
      ggtitle("Locus scores per epigenome") +
      coord_cartesian(ylim=c(-1,10))
    print(p)
  }
  
  dev.off()
  
  if (myargs$detailed) {
    pdf(file=paste0(myargs$output, ".loci.detailed.pdf"), width=myargs$width, height=myargs$height)
    plotsPerPage = 3
    numPages = ceiling(ncol(df) / plotsPerPage)
    for (i in 1:numPages) {
      starti = (i-1)*plotsPerPage + 1
      endi = min(i*plotsPerPage, ncol(df))
      plotlist = list()
      for (j in starti:endi) {
        subset.df = data.frame(eid=rownames(df), df[,j,drop=F])
        plot.df = tidyr::gather(subset.df, key=locus, value=score, 2:2)
        plot.df = plot.df[order(-plot.df$score),]
        plot.df$rank=1:nrow(plot.df)
        plot.df$topranks = paste0(1:nrow(plot.df), ". ", plot.df$eid)
        plot.df$topranks[7:nrow(plot.df)] = "other"
        plot.df$eid[!plot.df$eid %in% eids] = "other"
        plot.df$weight = 1
        plot.df$weight[plot.df$eid %in% eids] = 2
        p = ggplot(plot.df %>% dplyr::filter(eid == "other"), aes(x=locus, y=score, color=eid, size=weight, fill=topranks)) +
          geom_jitter(width=0.4) +
          geom_jitter(aes(x=locus, y=score, color=eid, size=weight),
                      data=plot.df %>% dplyr::filter(eid != "other"),
                      width=0.4) +
          scale_size(range=c(1,4)) + guides(size=FALSE) +
          theme_bw(9) + scale_color_discrete(name="EID") +
          theme(axis.text.x = element_text(angle=45, hjust=1)) +
          coord_cartesian(ylim=c(-2,9)) + theme(axis.text.y = element_blank())
        plotlist = append(plotlist, list(p))
      }
      grid.arrange(grobs=plotlist, ncol=plotsPerPage, top=textGrob("Locus scores per epigenome"))
    }
    dev.off()
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


main()

