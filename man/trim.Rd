\name{trim}
\alias{trim}
\title{Trim a word freqeuncy matrix}
\usage{
  trim(wfm, min.count = 5, min.doc = 5, sample.docs = NULL,
    sample.words = NULL, verbose = TRUE)
}
\arguments{
  \item{wfm}{a word count matrix with docs as rows and
  words as columns}

  \item{min.count}{minimum word frequency over all
  documents}

  \item{min.doc}{minimum number of documents containing at
  least one instance of a word}

  \item{sample.docs}{number of randomly chosen remaining
  word types to retain in the matrix.  If NULL, all
  documents are kept}

  \item{sample.words}{number of randomly chosen documents
  to retain in the matrix. If NULL, all words are kept}

  \item{verbose}{whether to generate a running commentary}
}
\value{
  a trimmed word frequency matrix
}
\description{
  Trim a word frequency matrix
}
\details{
  Removes words that do not both occur at least min.count
  times in at least min.doc documents, and optionally
  subsamples the remainder.
}
\author{
  Will Lowe
}

