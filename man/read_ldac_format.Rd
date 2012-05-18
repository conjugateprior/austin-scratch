\name{read_ldac_format}
\alias{read_ldac_format}
\title{Read an LDA-C file of word counts without word or document labels}
\usage{
  read_ldac_format(ldafile)
}
\arguments{
  \item{ldafile}{the file}
}
\value{
  a sparse matrix of word count information without word or
  documents
}
\description{
  Read LDA-C format file
}
\details{
  This functions reads a single file in LDA-C format.
}
\author{
  Will Lowe
}

