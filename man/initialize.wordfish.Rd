\name{initialize.wordfish}
\alias{initialize.wordfish}
\title{Initialize wordfish parameters}
\usage{
  initialize.wordfish(wfm)
}
\arguments{
  \item{wfm}{word count matrix with documents as rows and
  words as columns}
}
\value{
  a set of wordfish parameters: alpha, theta, psi and beta
}
\description{
  Initialize wordfish parameters
}
\details{
  This is mostly the Slapin and Proksch initialization
  code, which is in turn mostly the NOMINATE code involving
  an SVD of a double-centred word count matrix.  Use
  loglinear.initialize.wordfish instead.  It'll be much
  quicker.
}
\author{
  Will Lowe (from code by J. Slapin and S.-O. Proksch)
}

