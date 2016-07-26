\name{plot}
\docType{methods}
\alias{plot,OutlierDM-method}
\alias{plot,OutlierDM,missing-method}
\title{a plot-method for a "OutlierDM" object}
\description{
    This function provides a MA scatter plot with quantile regression based boxlplot. 
}
\usage{
    \S4method{plot}{OutlierDM}(x, y = NA, pch = 20, cex = 0.5, legend.use = TRUE, main, ...)
}
\arguments{
  \item{x}{
    fitted model object of class \code{\link{odm}}.  
  }
  \item{y}{
    the "y" argument is not used in the plot-method for "OutlierDM" object.
  }
  \item{pch}{
    a vector of plotting characters or symbols: see \code{\link[graphics:plot.default]{plot.defalut}}.  
  }
  \item{cex}{
    See \code{\link[graphics:plot.default]{plot.defalut}}.  
  }
  \item{legend.use}{  
    logical option for using legend box
  }
  \item{main}{
    main title for the plot
  }
  \item{...}{ \code{\link[graphics:plot.default]{plot.default}} arguments}
}
\details{
  This function is a method for the generic function \code{plot} for the S4 class 
  \code{OutlierDM}.  It can be invoked by calling print for an object of the 
  appropriate class, or directly by calling \code{plot.OutlierDM} regardless of 
  the class of the object. 
}
\seealso{
    \code{\link{odm}}
}
\keyword{subfunction}