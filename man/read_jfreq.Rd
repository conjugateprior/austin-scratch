\name{read_jfreq}
\alias{read_jfreq}
\title{Read a JFreq output folder}
\usage{
  read_jfreq(folder)
}
\arguments{
  \item{folder}{folder}
}
\value{
  a sparse Matrix object representing word frequencies with
  documents as rows and words as columns
}
\description{
  Read a JFreq output folder
}
\details{
  Reads the word frequency information from a folder
  created by JFreq.  All formats are supported.
}
\author{
  Will Lowe
}

