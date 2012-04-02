##' Read a JFreq output folder
##'
##' Reads the word frequency information from a folder created by JFreq.  All formats
##' are supported.
##' 
##' @title Read a JFreq output folder 
##' @param folder folder 
##' @export
##' @return a sparse Matrix object representing word frequencies with documents as rows and words as columns
##' @author Will Lowe
read_jfreq <- function(folder){
  csv <- grep("data.csv", dir(folder), value=TRUE)
  if (length(csv) > 0)
    m <- read.csv(file.path(folder, csv), encoding='UTF-8', row.names=1, header=TRUE)
  mtx <- grep("data.mtx", dir(folder), value=TRUE)
  if (length(mtx)>0)
    m <- read_mtx(folder, jfreq=TRUE)
  ldac <- grep("data.ldac", dir(folder), value=TRUE)
  if (length(ldac)>0)
    m <- read_ldac(folder, jfreq=TRUE)
  return(m)
}

##' Read a Matrix Market format file or folder
##'
##' If jfreq=TRUE read a JFreq output folder f in Matrix Market format, although
##' you could just use read_jfreq for that.  If jfreq=FALSE read a Matrix Market
##' file f with no row or column labels.  Note that jfreq determines whether f is
##' considered to be a folder or a file.  This is a thin wrapper around the Matrix
##' package's readMM function.
##' 
##' @title Read a Matrix Market format file or folder
##' @param f file or folder
##' @param jfreq whether to assume f is JFreq output
##' @return a sparse Matrix object representing word frequencies
##' @export
##' @author Will Lowe
read_mtx <- function(f, jfreq=TRUE){
  if (!jfreq)
    m <- read_mtx_format(f) ## just one file
  else {
    datafilename <- file.path(f, grep("data.mtx", dir(f), value=TRUE)) ## covers zipped ones
    m <- read_ldac_format(datafilename)
    try({rownames(m) <- read.csv(file.path(f, 'docs.csv'),
                                 fileEncoding='UTF8', header=FALSE)[,1]})
    try({colnames(m) <- read.csv(file.path(f, 'words.csv'),
                                 fileEncoding='UTF8', header=FALSE)[,1]})    
  }
  return(m)
}

##' Read a Matrix Market format file
##'
##' A thin wrapper around the Matrix package's readMM function.
##'
##' @title Read a Matrix Market format file
##' @param mtxfile a filename
##' @return a sparse Matrix object representing word frequencies
##' @export
##' @author Will Lowe
read_mtx_format <- function(mtxfile){
  readMM(mtxfile)
}

##' Read an LDA-C format file or folder
##'
##' When jfreq=TRUE read a JFreq output folder f in LDA-C format, although
##' you could just use read_jfreq for that.  When jfreq=FALSE read an LDA-C
##' file f with no row or column labels.  Note that jfreq determines whether f is
##' considered to be a folder or a file.  This function is a wrapper around read_ldac_format
##'
##' The LDA-C format is consists of a file in which line i is of the form
##' N [A:X]+ representing that the word with index A occurred X times in the i-th document.
##' Zero word counts are not represented so X>0 but A indexes into a word list starting
##' from 0.  N is the number of non-zero count words in the document and therefore
##' the number of colon-separated pairs in line i.
##' See the LDA-C documentation for details.
##'
##' Note: This implementation differs slightly from the lda package's read functions which look
##' for a data file and an accompanying word file for which A is the index.  It also does
##' not create a sparse matrix representation of the word count information.
##' 
##' @title Read an LDA-C format file or a folder of JFreq folder in LDA-C format
##' @param f file or folder
##' @param jfreq whether to assume f is JFreq output
##' @return a sparse Matrix object representing word frequencies
##' @export
##' @author Will Lowe
read_ldac <- function(f, jfreq=TRUE){
  if (!jfreq)
    m <- read_ldac_format(f) ## just one file
  else {
    ## a folder full of files
    datafilename <- file.path(f, grep("data.ldac", dir(f), value=TRUE)) ## covers zipped ones
    m <- read_ldac_format(datafilename)
    try({rownames(m) <- read.csv(file.path(f, 'docs.csv'),
                                 fileEncoding='UTF8', header=FALSE)[,1]})
    try({colnames(m) <- read.csv(file.path(f, 'words.csv'),
                                 fileEncoding='UTF8', header=FALSE)[,1]})    
  }
  return(m)
}

##' Read LDA-C format file
##'
##' This functions reads a single file in LDA-C format.
##' @title Read an LDA-C file of word counts without word or document labels
##' @param ldafile the file
##' @return a sparse matrix of word count information without word or documents
##' @export
##' @author Will Lowe
read_ldac_format <- function(ldafile){
  one <- scan(ldafile, what = "", sep = "\n")
  rs <- length(one) 
  two <- chartr(":", " ", one)
  three <- strsplit(two, " ", fixed = TRUE)
  ijx.len <- sum(unlist(lapply(three, function(x){ (length(x)-1)/2 })))
  i <- rep(0, ijx.len)
  j <- rep(0, ijx.len)
  x <- rep(0, ijx.len)
  start <- 1
  for (line in 1:rs){
    ll <- as.integer(three[[line]][-1])
    len.ll <- length(ll)
    end <- start + len.ll/2 - 1
    j[start:end] <- ll[seq(1, len.ll, by=2)] +1 ## lda-c format indexes from 0
    x[start:end] <- ll[seq(2, len.ll, by=2)]
    i[start:end] <- rep(line, len.ll/2)
    start <- end + 1
  }
  sparseMatrix(i, j, x=x)
}
##' Write a file for BMR/BXR
##'
##' This function writes out a list of word counts in sparse matrix format and a
##' dependent variable y in a form suitable for training the BMR/BXR family of
##' document classifiers.  The format is essentially the same as SVMLite.  Line i
##' is of the form Y [A:X]+ which represents that the i-th document is in class
##' number Y and that the word A occurs X times within it.  
##' 
##' @title Write a file for BMR/BXR
##' @param wfm a word frequency matrix with documents as rows and words as columns
##' @param y a numerical dependent variable, usually starting from 0
##' @param file name of the output file
##' @return nothing.  Used for the file saving side-effect
##' @export
##' @author Will Lowe
write_bmr <- function(wfm, y, file="data.bmr"){
  sink(file)
  for (i in 1:nrow(wfm)){
    nzero <- which(wfm[i,]>0) ## we index from zero, although this is not necessary
    val <- wfm[i,nzero]
    line <- paste(y[i], paste((nzero-1), ":", wfm[i,nzero], sep="", collapse=" "), "\n")
    cat(line)
  }
  sink()
}
##' Write a matrix into LDA-C format
##'
##' This function writes a sparse word frequency matrix into LDA-C format.  This involves
##' three files, called by default data.ldac, docs.csv and words.csv, and which contain
##' the word frequency information in LDA-C format, the document names, and the words
##' respectively.  See read_ldac for the details of the LDA-C format.
##' 
##' @title Write a matrix into LDA-C format
##' @param wfm a word frequency matrix with documents as rows and words as columns
##' @param folder the name of the folder containing the files. This will be created if necessary.
##' @param names filenames for the data, documents and word files
##' @return nothing. used for the file saving side-effect 
##' @export
##' @author Will Lowe
write_ldac <- function(wfm, folder, names=c("data.ldac", "docs.csv", "words.csv")){
  if (!file.exists(folder))
    dir.create(folder)
  ds <- file.path(folder, names[2])
  file.create(ds)
  ws <- file.path(folder, names[3])
  file.create(ws)  
  writeLines(rownames(wfm), ds)
  writeLines(colnames(wfm), ws)
  sink(file.path(folder, names[1]))
  for (i in 1:nrow(wfm)){
    nzero <- which(wfm[i,]>0) ## but index from zero
    val <- wfm[i,nzero]
    line <- paste(length(nzero), paste((nzero-1), ":", wfm[i,nzero], sep="", collapse=" "), "\n")
    cat(line)
  }
  sink()
}

##' Trim a word frequency matrix
##'
##' Removes words that do not both occur at least min.count times in total across at least min.doc documents,
##' and optionally subsamples the remainder.
##' 
##' @title Trim a word frequency matrix
##' @param wfm a word count matrix with docs as rows and words as columns
##' @param min.count minimum number of times a word occurs over all documents
##' @param min.doc minimum number of documents that must contain a word
##' @param sample.docs number of randomly chosen remaining word types to retain in the matrix.  If NULL, all documents are kept
##' @param sample.words number of randomly chosen documents to retain in the matrix. If NULL, all words are kept
##' @param verbose whether to generate a running commentary
##' @return a trimmed word frequency matrix
##' @export
##' @author Will Lowe
trim <- function(wfm, min.count=5, min.doc=5, sample.docs=NULL, sample.words=NULL, verbose=TRUE){
  
  N = nrow(wfm)
  V = ncol(wfm)
  
  rs1 <- which(colSums(wfm) >= min.count)
  if (verbose)
    cat("Words appearing less than", min.count, "times:", (V - length(rs1)), "\n")
  
  rs2 <- which(apply(wfm, 2, function(x){ sum(x>0) }) >= min.doc)
  if (verbose)
    cat("Words appearing in fewer than", min.doc, "documents:", (V - length(rs2)), "\n")
  
  tokeep <- intersect(rs1, rs2)
	
  if (length(tokeep)==0)
    stop("No words left after trimming!")
  
  if (!is.null(sample.words))
    tokeep.words <- sort(sample(tokeep, min(length(tokeep), sample.words)))
  else
    tokeep.words <- sort(tokeep) ## alphabetise
  
  if (!is.null(sample.docs))
    tokeep.docs <- sample(1:N, min(N, sample.docs))
  else
    tokeep.docs <- 1:N
  
  return(wfm[tokeep.docs, tokeep.words, drop=FALSE])
}


######## transformation functions

double.center <- function(X){
	rm <- rowMeans(X) 
	cm <- colMeans(X)
	X - outer(rm, cm, "+") + mean(X)
}


