wfm <- function(mat, word.margin=1){
    if (class(mat)=='character'){
        if (!file.exists(mat))
            stop(paste("The file", mat, "does not exist"))
        mat <- read.csv(file=file, row.names=1)
        if (word.margin==2)
            mat <- t(mat)
    }

    if (is.null(rownames(mat)) || is.null(colnames(mat))){
        stop("Cannot convert this object to a wfm: It must have row and column names")
    }

    ll <- as.matrix(mat)
    if (word.margin==1){
        dimnames(ll) <- list(words=rownames(mat), docs=colnames(mat))
    } else if (word.margin==2){
        dimnames(ll) <- list(words=colnames(mat), docs=rownames(mat))
    } else {
        stop("word.margin must be 1 (rows) or 2 (columns)")
    }

    return(ll)
}

is.wfm <- function(x){
    dnn <- names(dimnames(x))
    return(('words' %in% dnn) && ('docs' %in% dnn))
}

wordmargin <- function(x){
    if (!is.wfm(x))
        stop("Function not applicable to this object")
    wm <- ifelse(names(dimnames(x))[1]=='words', 1, 2)

    return(wm)
}

'wordmargin<-' <- function(x, value){
    if (!is.wfm(x))
        stop("Function not applicable to this object")
    if (value != 1 && value != 2)
        stop("New word margin must be 1 (rows) or 2 (columns)")

    if (wordmargin(x) != value){
        ## relabel the dimensions but don't change their contents
        if (value==1)
            dimnames(x) <- list(words=rownames(x), docs=colnames(x))
        else
            dimnames(x) <- list(docs=rownames(x), words=colnames(x))
    }

    return(x)
}

as.wfm <- function(mat, word.margin=1){
    ## rather speculative conversion from tm format here
    if (is(mat, 'TermDocumentMatrix')){
        return(wfm(as.matrix(mat), word.margin=1))
    }
    if (is(mat, 'DocumentTermMatrix')){
        return(wfm(as.matrix(mat), word.margin=2))
    }

    ## if it's a data.frame or a matrix this should work
    ## provided there are row and column names
    if (is.null(rownames(mat)) || is.null(colnames(mat)))
        stop("Cannot convert this object to a wfm.  It must have row and column names")

    ww <- wfm(mat, word.margin=word.margin)

    return(ww)
}

as.docword <- function(wfm){
    if (!is.wfm(wfm))
        stop("Function not applicable to this object")

    if (wordmargin(wfm)==1)
        return(t(wfm))
    else
        return(wfm)
}

as.worddoc <- function(wfm){
    if (!is.wfm(wfm))
        stop("Function not applicable to this object")

    if (wordmargin(wfm)==1)
        return(wfm)
    else
        return(t(wfm))
}

words <- function(wfm){
    if (wordmargin(wfm)==1)
        return(rownames(wfm))
    else
        return(colnames(wfm))
}

'words<-' <- function(x, value){
    if (length(words(x)) != length(value))
        stop("Replacement values are not the same length as the originals")

    if (wordmargin(x)==1)
        rownames(x) <- value
    else
        colnames(x) <- value
    return(x)
}

docs <- function(wfm){
    if (wordmargin(wfm)==1)
        return(colnames(wfm))
    else
        return(rownames(wfm))
}

'docs<-' <- function(x, value){
    if (length(docs(x)) != length(value))
        stop("Replacement value is not the same length as original")

    if (wordmargin(x)==1)
        colnames(x) <- value
    else
        rownames(x) <- value
    return(x)
}

trim <- function(wfm, min.count=5, min.doc=5, sample=NULL){
    ## eject words occurring less than min.count times and
    ## in fewer than min.doc documents
    ## then optionally sub sample words

    if (is.wfm(wfm))
        mY <- as.worddoc(wfm)
    else
        stop("Function not applicable to this object")

    rs1 <- which(rowSums(mY) < min.count)
    cat("Words appearing less than", min.count, "times:", length(rs1), "\n")

    rs2 <- which(apply(mY, 1, function(x){ (length(which(x>0)) < min.doc) }))
    cat("Words appearing in fewer than", min.doc, "documents:", length(rs2), "\n")

    togo <- union(rs1, rs2)
    res <- mY[-togo,]

    N <- NROW(res)
    if (!is.null(sample))
        voc <- sample(N, min(N,sample))
    else
        voc <- 1:N

    obj <- wfm(res[voc,], word.margin=1)

    return(obj)
}

wfm2lda <- function(wfm, dir=NULL, names=c("mult.dat", "vocab.dat")){
    m <- as.worddoc(wfm)
    v <- words(m)

    d <- list()
    for (i in 1:length(docs(m))){
        nzero <- which(m[,i]>0)
        d[[i]] <- t(matrix(as.integer(c(nzero-1, m[nzero, i])), ncol=2))
    }
    if (!is.null(dir))
        return(list(vocab=v, data=d))

    if (!file.exists(dir))
        stop(paste("Folder", dir, "does not exist"))
    lines <- rep("", length(d))
    for (i in 1:length(d)){
        lines[i] <- paste(NCOL(d[[i]]),
                          paste(d[[i]][1,], ":", d[[i]][2,], sep="", collapse=" "))
    }

    writeLines(lines, file.path(dir, names[1]))
    writeLines(v, file.path(dir, names[2]))
}

wfm2bmr <- function(y, wfm, filename){
    if (!is.null(y)){
        if (!(0 %in% y))
            stop("Dependent variable must index from 0")
        y <- as.numeric(y)+1 # to make 1 and 2 not 0 and 1
    }
    x <- as.worddoc(wfm)

    lines <- rep("", NCOL(x))
    for (i in 1:NCOL(x)){
        nonz <- which(x[,i] != 0.0)
        nonz.vals <- as.numeric(x[nonz,i])
        yy <- ifelse(!is.null(y), y[i], NULL)

        ldat <- paste(nonz, ":", nonz.vals, collapse=" ", sep="")
        lines[i] <- ifelse(is.null(y), ldat, paste(y[i], ldat))
    }

    writeLines(lines, filename)
}
