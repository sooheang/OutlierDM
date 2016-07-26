\name{toy}
\alias{toy}
\docType{data}
\title{
 an artificial dataset
}
\description{
  Artificial dataset from the simulation study of Eo et al. (2012). 
}
\usage{data("toy")}
\format{
  A data frame with 200 peptides (rows) and 15 samples (columns).
}
\source{
Eo, S-H, Pak, D, Choi, J, and Cho, H. (2012). Outlier Detection for Multiplicative High-throughput Data. \emph{BMC Research Notes}, \bold{5}, 1--6.
}
\examples{
data(toy)
## str(toy)
## pairs(log2(toy), pch = 20, cex = .7)
}
\keyword{datasets}
