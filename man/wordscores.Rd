\name{wordscores}
\alias{wordscores}
\title{Fit a wordscores model}
\usage{
  wordscores(wfm, scores, fit.fun = "classic.wordscores",
    verbose = FALSE, ...)
}
\arguments{
  \item{wfm}{a word count matrix, with documents as rows
  and columns as words}

  \item{scores}{reference scores for all documents}

  \item{fit.fun}{engine to use to fit the model}

  \item{verbose}{whether to give a running commentary}

  \item{...}{further arguments passed to engine specified
  by fit.fun}
}
\value{
  a wordscores model
}
\description{
  Fit a wordscores model
}
\details{
  This function computes wordscores using the algorithm
  described in Laver and Garry (2003) and implemented as
  classic.wordscores. To score 'virgin', i.e out-of-sample
  documents, use the predict function. Words that do not
  occur in the document set are quietly removed before
  wordscores are computed.
}
\author{
  Will Lowe
}

