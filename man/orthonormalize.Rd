\name{orthonormalize}
\alias{orthonormalize}
\title{Orthonormalise the columns of a matrix}
\usage{
  orthonormalize(u, basis = FALSE, norm = FALSE)
}
\arguments{
  \item{u}{a matrix}

  \item{basis}{whether to create a basis too}

  \item{norm}{whether to orthonormalise too}
}
\value{
  a matrix the same size with orthogonal columns
}
\description{
  Orthonormalize the columns of a matrix
}
\details{
  This function orthogonalises the columns of a matrix u
  using the Gram-Schmidt method.  It is essentially the
  orthonormalisation function from the far package by
  Julien Damon.
}
\author{
  Will Lowe (original code by Damon Julien Guillas Serge)
}

