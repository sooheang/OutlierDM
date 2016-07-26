\name{odm.control}
\alias{odm.control}
\title{Control tuning parameters for "OutlierDM" object}
\description{
	various parameters that control aspects of the "OutlierDM" object
}
\usage{
     odm.control(pair.cre = 1, dist.mthd = "median", 
     	Lower = .25, Upper = .75, trans = "log2", 
     	centering = TRUE, projection.type = "PCA", lbda = 1, 
		nonlin.method = "L-BFGS-B", nonlin.SS = "AsymOff", 
		nonlin.Frank = c(2, -8, 0, 1), ncl = 2, cri.pval = 0.05 )
}
\arguments{
   \item{pair.cre}{a scalar parameter to specify the minimum number of pairs, used in "type = pair".}
   \item{dist.mthd}{a distance parameter used in "type = diff". you can choose one of "median", "mean" or and so on.}
   \item{Lower}{a criterion for lower quantile value used to construct boxplot}
   \item{Upper}{a criterion for upper quantile value used to construct boxplot}
   \item{trans}{a parameter for a logarithm and exponential transformation. If a log 2 transformation is needed, set "trans = log2". If no transformation is needed, set "trans = FALSE". }
   \item{centering}{a logical parameter for the status of centering. If "centering = TRUE", data are centered by its column means.}
   \item{projection.type}{a parameter to determine a type of projection methods. Choose one of "naive", "pca", and "robust".}
   \item{lbda}{a criterion about lambda used for nonlinear quantile regression.}
   \item{nonlin.method}{a parameter to determine a type of methods used for nonlinear quantile regression. choose one of "L-BFGS-B" and "BFGS". Default is "L-BFGS-B".}
   \item{nonlin.SS}{a parameter to determine a type of structure used for nonlinear quantile regression. choose one of "Frank", "Self", "Asym" and "AsymOff". Default is "AsymOff", Asymptotic Regression Model with an Offset.}
   \item{nonlin.Frank}{a structure parameter used for Frank copula model. Gain c(df, delta, mu, sigma) in the Frank copula formula}
   \item{ncl}{A parameter to determine the number of cores used in parallel computing. A default value is 2.}
   \item{cri.pval}{A criteria for the p-value of the Grubbs test. A default value is 0.05.}
}
\seealso{
    \code{\link{odm}}
}
\keyword{odm}


