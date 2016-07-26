\name{rgrubbs.test}
\alias{rgrubbs.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recursive Grubbs test for an outlier detection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function works to detect outlying observation given one peptide by using the Grubbs test recursively.
}
\usage{
rgrubbs.test(x, cri = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a peptide}
  \item{cri}{a criteria for the p-value}
}
\details{
It is a recursive version of the Grubbs test to detect outlying observations assuming that a peptide is given.
}
\references{
%% ~put references to the literature/web site here ~
Cho, H and Eo, S-H. (2015). Outlier Detection for Mass Spectrometry Data, \emph{Submitted}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{odm}}
}
\examples{
  data(lcms3)
  rgrubbs.test(log2(lcms3[100,]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subfunction }
