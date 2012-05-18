\name{summary.insample.wordscores}
\alias{summary.insample.wordscores}
\title{Summarise an 'in sample 'wordscores model}
\usage{
  \method{summary}{wordscores} (object, ...)
}
\arguments{
  \item{object}{insample.wordscores model}

  \item{...}{Not used}
}
\value{
  A data.frame with document statistics
}
\description{
  Summarise an 'in sample' wordscores model
}
\details{
  Generates some minimal statistics about the document set
  from which wordscores have been calculated, and the
  estimated positions of unscored documents
}
\author{
  Will Lowe
}

