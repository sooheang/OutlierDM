\name{OutlierDM-class}
\Rdversion{1.1}
\docType{class}
\alias{OutlierDM-class}
\alias{input}
\alias{output}
\alias{outliers}
\alias{show,OutlierDM-method}
\alias{input,OutlierDM-method}
\alias{output,OutlierDM-method}
\alias{outliers,OutlierDM-method}
\title{Class \code{"OutlierDM"}}
\description{
	A \code{S4} class for the \code{OutlierDM} package
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("OutlierDM", ...)}. See following information about slots.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{evaluated function call}
    \item{\code{raw.data}:}{data to be used in the fitted model}
    \item{\code{res}:}{a data.frame including the information about the fitted model. It consists of several columns including \code{outlier, M, A, Q3, Q1, UB} and \code{LB}.}
    \item{\code{x.pair}:}{a list including the information of the pairwise outlierD algorithm}
    \item{\code{k}:}{a scalar parameter for constructing boxplot used in the fitted models}
    \item{\code{outlier}:}{a boolean matrix for outlier information}
    \item{\code{n.outliers}:}{a scalar value that denotes the number of outliers to be detected by the fitted model.}
    \item{\code{method}:}{the type of method used in the fitted model}
    \item{\code{type}:}{the type of algorithm used in the fitted model}
    \item{\code{contrl.para}:}{a list including information about tuning parameters}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "OutlierDM")}: Same as the show method without the optional arguments}
    \item{plot}{\code{signature(x = "OutlierDM", y = "missing")}: Plot an object.}
    \item{input}{\code{signature{object = "OutlierDM"}}: Show an input data set}
    \item{output}{\code{signature{object = "OutlierDM"}}: Show the result}
    \item{outliers}{\code{signature{object = "OutlierDM"}}: Show the candidate outliers}
	 }
}
\examples{
  showClass("OutlierDM")
}
\seealso{
  \code{\link{odm}}
}
\keyword{odm}
