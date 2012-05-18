\name{classic.wordscores}
\alias{classic.wordscores}
\title{Fit a wordscores model}
\usage{
  classic.wordscores(wfm, scores, verbose)
}
\arguments{
  \item{wfm}{a word count matrix with documents as rows and
  words as columns}

  \item{scores}{reference document scores}

  \item{verbose}{ignored}
}
\value{
  a fitted wordscores model
}
\description{
  Fit a Wordscores model
}
\details{
  Fit a Wordscores model using the algorithm in Laver and
  Garry. This means tossing out the words that never occur
  and constructing wordscores. Consider using
  correspondence analysis instead!
}
\author{
  Will Lowe
}

