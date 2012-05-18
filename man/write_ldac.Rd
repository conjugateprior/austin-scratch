\name{write_ldac}
\alias{write_ldac}
\title{Write a matrix into LDA-C format}
\usage{
  write_ldac(wfm, folder,
    names = c("data.ldac", "docs.csv", "words.csv"))
}
\arguments{
  \item{wfm}{a word frequency matrix with documents as rows
  and words as columns}

  \item{folder}{the name of the folder containing the
  files. This will be created if necessary.}

  \item{names}{filenames for the data, documents and word
  files}
}
\value{
  nothing. used for the file saving side-effect
}
\description{
  Write a matrix into LDA-C format
}
\details{
  This function writes a sparse word frequency matrix into
  LDA-C format.  This involves three files, called by
  default data.ldac, docs.csv and words.csv, and which
  contain the word frequency information in LDA-C format,
  the document names, and the words respectively.  See
  read_ldac for the details of the LDA-C format.
}
\author{
  Will Lowe
}

