\name{predict.wordfish}
\alias{predict.wordfish}
\title{Predictions for new documents.}
\usage{
  \method{predict}{wordfish} (object, newdata = NULL,
    se.fit = FALSE, interval = c("none", "confidence"),
    level = 0.95, ...)
}
\arguments{
  \item{object}{fitted wordfish model}

  \item{newdata}{new data or the original data used to fit
  the model if NULL}

  \item{se.fit}{whether standard errors are reported}

  \item{interval}{whether confidence intervals are
  reported}

  \item{level}{coverage as a proportion}

  \item{...}{Not used}
}
\value{
  a data.frame containing predicted positions for each
  document in newdata
}
\description{
  Predictions for new documents
}
\details{
  This function provides predicted values and uncertainty
  estimates for new documents assuming that other model
  coefficients are known.  The multinomial formulation is
  used to create analytic measures of uncertainty.  A
  better version of this function would take into account
  uncertainty in the word parameter values when
  constructing the latter.
}
\author{
  Will Lowe
}

