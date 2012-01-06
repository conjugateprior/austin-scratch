\name{print.wordfish}
\alias{print.wordfish}
\title{Show a wordfish model}
\usage{
  \method{print}{wordfish} (x,
    digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{wordfish model}

  \item{digits}{how many decimal places to use}

  \item{...}{Not used.}
}
\value{
  Nothing. Used for printing side-effect
}
\description{
  Show a wordfish model
}
\details{
  Shows the fitted document positions and intervals for a
  fitted wordfish model.
}
\author{
  Will Lowe
}

