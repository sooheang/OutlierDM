\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
\name{NEWS}
\title{OutlierDM News}
\encoding{UTF-8}

% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.0 2014-12-31}{
  \itemize{
  	\item Add four classical outlier detection approaches for large samples
    \itemize{
      \item Zscore: Z-score based criterion
      \item grubbs: recurive Grubbs test
      \item iqr: interquartile range (IQR) criterion
      \item siqr: semi-interquartile range (SIQR) criterion
    }
    \item Add new plotting function called oneplot() designed to detect outliers given a peptide
    \item Add projection pursuit PCA for the projection method
    \item Add toy dataset consisting with 200 peptides and 15 samples
    \item Bug fixed on manuals
    \item Redesign the output of the fitted object
    \item ToDo list
    \itemize{
      \item Connect to a web-based engine
      \item Redesign plotting functions
      \item Expand the algorithm to Empirical Bayesian approaches
    }
  }
}

% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0-0 2014-03-22}{
  \itemize{
    \item Publish an official version on CRAN
    \item ToDo 1: Add new routine for finding optimal K
    \item ToDo 2: Add new algorithm and routine using Empirical Bayesian approach
  }
}
