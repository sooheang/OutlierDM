\name{OutlierDM-package}
\alias{OutlierDM}
\alias{OutlierDM-package}
\docType{package}
\title{
	Functions for detecting outlying observations in (multiple) replicated high-throughput data
}
\description{
  This package provides outlier detection algorithms for multiple replicated high-throughput data, especially in the field of mass spectrometry.
}
\details{
\tabular{ll}{
Package: \tab OutlierDM\cr
Type: \tab Package\cr
Version: \tab 1.1-0\cr
Date: \tab 2014-12-31\cr
License: \tab GPL version 3\cr
LazyLoad: \tab no\cr
}
}
\author{
Soo-Heang Eo <eo.sooheang@gmail.com> and HyungJun Cho <hj4cho@korea.ac.kr>

Maintainer: Soo-Heang Eo <eo.sooheang@gmail.com>
}
\references{
Cho H, and Eo, S-H. (2014). Outlier Detection for Mass-Spectrometry Data. \emph{Submitted}. 

Eo, S-H, Pak, D, Choi, J, and Cho, H. (2012). Outlier Detection for Multiplicative High-throughput Data. \emph{BMC Research Notes}, \bold{5}, 1--6.

Cho et al. (2008). OutlierD: an R package for outlier detection using quantile regression on mass spectrometry data. \emph{Bioinformatics}, \bold{24}(6), 882--884.

Min et al. (2007). Ultrahigh-pressure dual online solid phase extraction/capillary reverse-phase liquid chromatography/tandem mass spectrometry (DO-SPE/cRPLC/MS/MS): A versatile separation platform for high-throughput and highly sensitive proteomic analyses. \emph{ELECTROPHORESIS}, \bold{28}, 1012--1021. 
}

\examples{
\dontrun{
data(lcms3)

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
\seealso{
    \code{\link{odm}}, \code{\link{odm.control}}, \code{\link[quantreg:rq]{quantreg}}
}
\keyword{package}