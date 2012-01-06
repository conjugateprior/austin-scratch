\name{read.mtx.format}
\alias{read.mtx.format}
\title{Read a Matrix Market format file}
\usage{
  read.mtx.format(mtxfile)
}
\arguments{
  \item{mtxfile}{a filename}
}
\value{
  a sparse Matrix object representing word frequencies
}
\description{
  Read a Matrix Market format file
}
\details{
  A thin wrapper around the Matrix package's readMM
  function.
}
\author{
  Will Lowe
}

