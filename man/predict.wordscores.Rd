\name{predict.wordscores}
\alias{predict.wordscores}
\title{Score new document with a wordscores model}
\usage{
  \method{predict}{wordscores} (object, newdata = NULL,
    rescale = c("lbg", "none"), level = 0.95, ...)
}
\arguments{
  \item{object}{wordscores model}

  \item{newdata}{A new word coutn matrix.  If NULL, the
  original matrix}

  \item{rescale}{whether to apply the rescaling described
  in Laver Benoit and Garry 2003 which fixes the new
  documents' variance to be the same as the original
  documents.}

  \item{level}{Notional coverage as a proportion}

  \item{...}{Not used.}
}
\value{
  A data.frame containing predited new document positions
}
\description{
  Score new documents with a wordscores model
}
\details{
  This function uses the wordscores from a fitted model to
  estimate document positions for new documents, a.k.a.
  virgin documents.  If newdata=NULL the original
  'reference' documents are used (their scores are not in
  general recovered due to shrinkage and other factors...)

  Standard errors are computed as per the original paper,
  which may or may not be sensible; the idea of standard
  errors in the absence of an explicit probability model is
  unclear.
}
\author{
  Will Lowe
}

