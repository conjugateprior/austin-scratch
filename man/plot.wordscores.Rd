\name{plot.wordscores}
\alias{plot.wordscores}
\title{Plot wordscores}
\usage{
  \method{plot}{wordscores} (x, ...)
}
\arguments{
  \item{x}{a wordscores model}

  \item{...}{Graphics parameters passed on to dotchart}
}
\value{
  Nothing. Used for plotting side effect
}
\description{
  Plot the wordscores from a model
}
\details{
  Makes a dotchart of wordscores, order in size. Not
  especially informative, except for viewing the
  document-driven spikes described in Lowe 2008.
}
\author{
  Will Lowe
}

