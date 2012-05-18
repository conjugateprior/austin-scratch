\name{read_mtx}
\alias{read_mtx}
\title{Read a Matrix Market format file or folder}
\usage{
  read_mtx(f, jfreq = TRUE)
}
\arguments{
  \item{f}{file or folder}

  \item{jfreq}{whether to assume f is JFreq output}
}
\value{
  a sparse Matrix object representing word frequencies
}
\description{
  Read a Matrix Market format file or folder
}
\details{
  If jfreq=TRUE read a JFreq output folder f in Matrix
  Market format, although you could just use read_jfreq for
  that.  If jfreq=FALSE read a Matrix Market file f with no
  row or column labels.  Note that jfreq determines whether
  f is considered to be a folder or a file.  This is a thin
  wrapper around the Matrix package's readMM function.
}
\author{
  Will Lowe
}

