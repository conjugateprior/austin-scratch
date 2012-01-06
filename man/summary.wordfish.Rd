\name{summary.wordfish}
\alias{summary.wordfish}
\title{Summarise a wordfish model}
\usage{
  \method{summary}{wordfish} (object, level = 0.95, ...)
}
\arguments{
  \item{object}{fitted wordfish model}

  \item{level}{coverage as a proportion}

  \item{...}{Not used}
}
\value{
  A lightweight summary object containing the model call
  and position estimates
}
\description{
  Summarise a wordfish model
}
\details{
  Summarises the fitted document positions theta (row
  scores) with confidence intervals
}
\author{
  Will Lowe
}

