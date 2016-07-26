\name{oneplot}
\alias{oneplot, OutlierDM, numeric-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Draw a dot-plot for a selected peptide
}
\description{
This function draws a dot plot for a selected peptide based on the \code{OutlierDM} object.
}
\usage{
    \S4method{oneplot}{OutlierDM}(object,i=1,pick=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted object}
  \item{i}{ a row number in order or drawing a dow plot}
  \item{pick}{the number of locators to denote the names of the peptide}
}
\seealso{
\code{\link{odm}}
}
\examples{
  data(lcms3)
  fit = odm(lcms3, type = "grubbs")
  oneplot(fit, i = 100)

  # Add row name
  #oneplot(fit, i = 100, pick = 1)
}
\keyword{subfunction}
