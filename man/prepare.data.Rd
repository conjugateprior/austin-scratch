\name{prepare.data}
\alias{prepare.data}
\title{Tidy up count matrix before estimation}
\usage{
  prepare.data(wfm, verbose)
}
\arguments{
  \item{wfm}{word count matrix: document as rows, words as
  columns}

  \item{verbose}{whether to offer running commentary}
}
\value{
  a word count matrix
}
\description{
  Tidy up count matrix before estimation
}
\details{
  This function removes words that do not occur in the
  corpus.  Mostly useful for when bootstrap functions make
  degenerate word count matrices.
}
\author{
  Will Lowe
}

