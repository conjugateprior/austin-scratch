\name{plot.wordfish}
\alias{plot.wordfish}
\title{Plot a wordfish model}
\usage{
  \method{plot}{wordfish} (x, truevals = NULL,
    level = 0.95, pch = 20, ...)
}
\arguments{
  \item{x}{fitted wordfish model}

  \item{truevals}{true positions for comparison}

  \item{level}{coverage as a proportion}

  \item{pch}{plotting character}

  \item{...}{parameters passed to dotchart}
}
\value{
  Nothing. Used for plotting side effect
}
\description{
  Plot a wordfish model
}
\details{
  Plots document postions, associated measures of
  uncertainty, and optionally true positions if they are
  available.
}
\author{
  Will Lowe
}

