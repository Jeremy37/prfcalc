#!/nfs/users/nfs_j/js29/R-3.0.3/bin/Rscript
library(RobustRankAggreg)
source(ppaFunctions.R)

myargs <- NULL


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
  
  options(datatable.fread.datatable=T)
  options(datatable.verbose = F)
  options(datatable.showProgress = F)
  
  df <- fread(myargs$file, stringsAsFactors=T, header=T, sep='\t', na.strings=c("NA","N/A",""))
  df.t <- df
  
  if (!is.null(myargs$segmentCol)) {
    df.t$segment <- df[,myargs$segmentCol,with=F][[1]]
  }
  if (!is.null(myargs$pCol)) {
    df.t$Pval <- as.double(df[,myargs$pCol,with=F][[1]])
  }
  if (!is.null(myargs$F)) {
    df.t$F <- as.double(df[,myargs$F,with=F][[1]])
  }
  if (!is.null(myargs$N)) {
    df.t$N = as.integer(myargs$N)
  }
  if (!is.null(myargs$debug)) {
    myargs$debug = T
  } else {
    myargs$debug = F
  }
  
  # Check that required columns are present
  if (! ("Pval" %in% colnames(df.t) & "N" %in% colnames(df.t) & "F" %in% colnames(df.t))) {
    stop("Missing required column: 'Pval', 'N', and 'F' columns must be present.")
  }
  
  df[,"segment"] <- getSegmentNaivePPAs(df.t)
  
  if (myargs$verbose) {
    write.table(df, "", row.names=F, col.names=T, quote=F, sep="\t")
  }
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


