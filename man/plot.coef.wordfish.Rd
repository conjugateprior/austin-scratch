\name{plot.coef.wordfish}
\alias{plot.coef.wordfish}
\title{Plot model parameters}
\usage{
  \method{plot}{coef.wordfish} (x, pch = 20, psi = TRUE,
    ...)
}
\arguments{
  \item{x}{wordfish model parameters}

  \item{pch}{plotting character: defaults to fat dots}

  \item{psi}{whether to plot beta or beta against psi}

  \item{...}{Not used}
}
\value{
  Nothing. Used for plotting side-effect
}
\description{
  Plot model parameters
}
\details{
  Plots word parameter extracted using coef with
  form='poisson', in the familiar 'Eiffel Tower' plot (beta
  against psi) or a dotchart of beta if psi=FALSE.
}
\author{
  Will Lowe
}

