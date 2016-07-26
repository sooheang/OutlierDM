\name{odm}
\alias{odm}
\title{Outlier dectection for (multiple) replicated high-throughput data}
\description{
	This function provides some routines for detecting outlying observations (peptides) for multiple replicated high-throughput data, especially in mass-spectrometry by using quantile regression-based boxplot algorithms.
}
\usage{
  odm(x, k=1.5, method= c("linear", "nonlin", "constant", "nonpar"), 
  type = c("proj", "diff", "pair", "grubbs", "dixon", "iqr",  "siqr", "Zscore"), ...)
}

\arguments{
  \item{x}{
    data vectors or matrices. These can be given as named arguments. If the number of predictors is 2, x1 describes one n-by-1 vector for data and x2 describes the other n-by-1 vector for data (n= number of peptides, proteins, or genes)
  }
  \item{k}{
    non-negative tuning parameter for the outlier detection algorithm. For the sake of IQR-based algorithms such as \code{'iqr'}, \code{'proj'}, \code{'diff'}, and \code{'pair'}, it works in the formula of Q1-k*IQR and Q3+k*IQR, IQR=Q3-Q1. For \code{'Zscore'}, it works for the \code{'k'} in |Z| > k. A default value is 1.5
  }
  \item{method}{type of quantile regression methods used in an outlier detection algorithm. Use one of \code{'constant'}, \code{'linear'}, \code{'nonlin'}, and \code{'nonpar'} which mean constant, linear, non-linear, and non-parametric quantile regression in order. For more details, see the \code{quantreg} package.}
  \item{type}{type of criterion for detecting outlying observations. You can select one of \code{'Zscore'}, \code{'iqr'}, \code{'dixon'}, \code{'grubbs'}, \code{'pair'}, \code{'diff'}, and \code{'proj'} algorithms. \cr
    Zscore: Z-score based criterion \cr
    iqr: Interquartile range (IQR) criterion \cr
    siqr: Semi-interquartile range (IQR) criterion \cr
  	dixon: Dixon's Q-test (Dixon, 1950; 1951) \cr
  	grubbs: Grubbs test (Grubbs, 1950; 1969) \cr
  	pair: Pariwise OutlierD algorithm (Cho et al., 2008; Eo et al., 2012) \cr
  	diff: Difference approach for multiplicative high-throughput data (Eo et al., 2012)\cr
  	proj: Boxplot-based MA plot with projections (Eo et al., 2012) \cr	
  	}
  \item{...}{
    minor tuning parameters used in odm.control(). See \code{\link{odm.control}}.
  }
}
\value{
   \item{\code{call}:}{evaluated function call}
    \item{\code{raw.data}:}{data to be used in the fitted model}
    \item{\code{res}:}{a data.frame including the information about the fitted model. It consists of several columns including \code{outlier, M, A, Q3, Q1, UB} and \code{LB}.}
    \item{\code{x.pair}:}{Object of class \code{"list"} }
    \item{\code{k}:}{a scalar parameter for constructing boxplot used in the fitted models}
    \item{\code{n.outliers}:}{a scalar value that denotes the number of outliers to be detected by the fitted model.}
    \item{\code{method}:}{the type of method used in the fitted model}
    \item{\code{type}:}{the type of algorithm used in the fitted model}
    \item{\code{contrl.para}:}{a list including information about tuning parameters}
}
\seealso{
    \code{\link{OutlierDM-package}} to provide the general information about the OutlierDC package\cr
	\code{\link{OutlierDM-class}} to provide the information about the \code{"OutlierDM"} class \cr
    \code{\link{odm.control}} to control tuning parameters \cr
}
\references{

Cho H, Lee JW, Kim Y-J, et al. (2008) OutlierD: an R package for outlier detection using quantile regression on mass spectrometry data. Bioinformatics 24:882–884. 

Dixon WJ (1950) Analysis of extreme values. Ann Math Statistics 21:488–506. 

Dixon WJ (1951) Ratios involving extreme values. Ann Math Statistics 22:68–78. 

Grubbs FE (1950) Sample criteria for testing outlying observations. Ann Math Statistics 21:27–58. 

Grubbs FE (1969) Procedures for detecting outlying observations in samples. Technometrics 11:1–21. 

1. Eo S-H, Pak D, Choi J, Cho H (2012) Outlier detection using projection quantile regression for mass spectrometry data with low replication. BMC Res Notes. doi: 10.1186/1756-0500-5-236
}
\examples{
  \dontrun{

    ## Load datasets
    data(lcms3)
    
    ## Load 
    ## Fit projection approaches
    fit.proj.const <- odm(lcms3, method="constant")
    fit.proj.linear <- odm(lcms3, method="linear")
    fit.proj.nonlin <- odm(lcms3, method="nonlin")
    fit.proj.nonpara <- odm(lcms3, method="nonpar", lbda = 1)
    
    par(mfrow = c(2,2))
    plot(fit.proj.const, main = "Constant")
    plot(fit.proj.linear, main = "Linear")
    plot(fit.proj.nonlin, main = "NonLinear")
    plot(fit.proj.nonpara, main = "Nonparametric")

    ## Fit pairwise OutlierD algorithm
    fit0 <- odm(lcms3, type = "pair")
    plot(fit0)

  }
}
\keyword{odm}


