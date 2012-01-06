\name{sim.wordfish}
\alias{sim.wordfish}
\title{simulate wordfish data}
\usage{
  sim.wordfish(docs = 10, vocab = 20)
}
\arguments{
  \item{docs}{how many documents to create}

  \item{vocab}{how many words to create}
}
\value{
  a named list containing simulated data and parameter
  values
}
\description{
  Simulate wordfish data
}
\details{
  Simulates wordfish data where psi=2, alpha[1]=0,
  alpha[2:docs]=1, theta consists of docs evenly spaced
  values with mean 0 and variance 1, and beta consists
  vocab evenly spaced values between -1.5 and 1.5. A
  simulated data matrix Y and these parameter values are
  returned a named list.
}
\author{
  Will Lowe
}

