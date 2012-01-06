################### Package ################################

#' Austin does things with words
#'
#' Austin implements Wordfish, Wordscores, and related text scaling models.  
#' It also provides utility function for working with various sparse 
#' data formats used in text analysis, e.g. LDA-C, Matrix Market, and variants
#' on the SVMLite format.  
#'
#' @import Matrix numDeriv
#' @author Will Lowe \email{will.lowe@@uni-mannheim.de}
#' @docType package
#' @name austin 
#' @aliases austin package-austin
NULL

#' German party manifestos
#'
#' A random sample of words from German party manifestos from 1990-2005.
#'
#' @references J. Slapin and S.-O. Proksch (2008) 'A scaling model for estimating time-series party positions from texts' American Journal of Political Science 52(3), 705-722
#' @name demanif
#' @docType data
#' @author Will Lowe, from data collected by Slapin and Proksch
#' @keywords data
NULL

#' German party manifestos, economic issues
#'
#' A random sample of words from sections of German party manifestos concerning economic 
#' issues from 1990-2005.
#'
#' @references J. Slapin and S.-O. Proksch (2008) 'A scaling model for estimating time-series party positions from texts' American Journal of Political Science 52(3), 705-722
#' @name demanif.econ
#' @docType data
#' @author Will Lowe, from data collected by Slapin and Proksch
#' @keywords data
NULL

#' German party manifestos, social issues
#'
#' A random sample of words from sections of German party manifestos concerning social 
#' issues from 1990-2005.
#'
#' @references J. Slapin and S.-O. Proksch (2008) 'A scaling model for estimating time-series party positions from texts' American Journal of Political Science 52(3), 705-722
#' @name demanif.soc
#' @docType data
#' @author Will Lowe, from data collected by Slapin and Proksch
#' @keywords data
NULL

#' German party manifestos, foreign policy
#'
#' A random sample of words from sections of German party manifestos concerning foreign policy 
#' issues from 1990-2005.
#'
#' @references J. Slapin and S.-O. Proksch (2008) 'A scaling model for estimating time-series party positions from texts' American Journal of Political Science 52(3), 705-722
#' @name demanif.foreign
#' @docType data
#' @author Will Lowe, from data collected by Slapin and Proksch
#' @keywords data
NULL

#' Irish no-confidence debate
#'
#' Speeches from a no-confidence debate from Ireland [DATE]
#'
#' @name daildata
#' @docType data
#' @author Will Lowe, from data collected by Ken Benoit
#' @keywords data
NULL

#' Irish budget debate
#'
#' Speeches from the budget debate in Ireland 2009.
#'
#' @name iebudget2009
#' @docType data
#' @author Ken Benoit
#' @keywords data
NULL

#' Irish budget data
#'
#' Speeches from the budget debate in Ireland [DATE]
#'
#' @name budgetdata
#' @docType data
#' @author Ken Benoit
#' @keywords data
NULL

#' Interest groups and the European Commission
#'
#' Word counts from interest group statements and two versions of the 
#' European commission's proposals to reduce CO2 emmissions, from 2007.
#'
#' \code{comm1} and \code{comm2} are the commission's proposal before and after
#' the proposals of the interest groups respectively.
#'
#' @name interestgroups
#' @docType data
#' @author Heike Kluever
#' @references H. Kluever (2009) 'Measuring influence group influence using quantitative text analysis' European Union Politics 11:1.
#' @keywords data
NULL

#' UK party manifestos
#'
#' Word counts from UK party manifestos, in 1992 and 1997.
#'
#' @name ukmanif
#' @docType data
#' @author Michael Laver, Ken Benoit, and John Garry
#' @references Laver, M. Benoit, K. and Garry, J (2003) Extracting policy positions from political texts using words as data. American Political Science Review 97:2
#' @keywords data
NULL

#' Fake wordscores data
#'
#' Fake data for showing the wordscores algorithm.  Documents R1-R5 are 'reference' 
#' with scores: -1.5, -0.75, 0, 0.75, 1.5 and
#' V1 is 'virgin' with true score 0.45.
#'
#' @name lbg
#' @docType data
#' @author Michael Laver, Ken Benoit, and John Garry
#' @references Laver, M. Benoit, K. and Garry, J (2003) Extracting policy positions from political texts using words as data. American Political Science Review 97:2
#' @keywords data
NULL




