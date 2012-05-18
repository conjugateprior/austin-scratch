\name{read_ldac}
\alias{read_ldac}
\title{Read an LDA-C format file or a folder of JFreq folder in LDA-C format}
\usage{
  read_ldac(f, jfreq = TRUE)
}
\arguments{
  \item{f}{file or folder}

  \item{jfreq}{whether to assume f is JFreq output}
}
\value{
  a sparse Matrix object representing word frequencies
}
\description{
  Read an LDA-C format file or folder
}
\details{
  When jfreq=TRUE read a JFreq output folder f in LDA-C
  format, although you could just use read_jfreq for that.
  When jfreq=FALSE read an LDA-C file f with no row or
  column labels.  Note that jfreq determines whether f is
  considered to be a folder or a file.  This function is a
  wrapper around read_ldac_format

  The LDA-C format is consists of a file in which line i is
  of the form N [A:X]+ representing that the word with
  index A occurred X times in the i-th document. Zero word
  counts are not represented so X>0 but A indexes into a
  word list starting from 0.  N is the number of non-zero
  count words in the document and therefore the number of
  colon-separated pairs in line i. See the LDA-C
  documentation for details.

  Note: This implementation differs slightly from the lda
  package's read functions which look for a data file and
  an accompanying word file for which A is the index.  It
  also does not create a sparse matrix representation of
  the word count information.
}
\author{
  Will Lowe
}

