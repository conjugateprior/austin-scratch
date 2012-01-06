\name{write.bmr}
\alias{write.bmr}
\title{Write a file for BMR/BXR}
\usage{
  write.bmr(wfm, y, file = "data.bmr")
}
\arguments{
  \item{wfm}{a word frequency matrix with documents as rows
  and words as columns}

  \item{y}{a numerical dependent variable, usually starting
  from 0}

  \item{file}{name of the output file}
}
\value{
  nothing.  Used for the file saving side-effect
}
\description{
  Write a file for BMR/BXR
}
\details{
  This function writes out a list of word counts in sparse
  matrix format and a dependent variable y in a form
  suitable for training the BMR/BXR family of document
  classifiers.  The format is essentially the same as
  SVMLite.  Line i is of the form Y [A:X]+ which represents
  that the i-th document is in class number Y and that the
  word A occurs X times within it.
}
\author{
  Will Lowe
}

