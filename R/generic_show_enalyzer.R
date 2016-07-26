#' Show an Enalyzer Object
#'
#' This function displays the \code{Enalyzer} object by printing by callging 
#'
#'@export
#'@import tidyr
#'@importFrom lubridate
#' @param formula symbolic description of the model of type y ~ x | z
#' @param data argument controlling fomrmula processing via \link{model.freame}
#' @param subset argument controlling fomrmula processing via \link{model.freame}
#' @param na.action argument controlling fomrmula processing via \link{model.freame}
#' @param weights optional numeric vector of case weights
#' @param offset optional numericvector with an a priori known component 
#' @param type type for the analysis (under construction)
#' @param control a list of minor arguments spcified via \link{enalyzerControl}.
#' @param verbose tracing the procedure 
#' @param ... list artuments (under construction)

setGeneric("show")
setMethod("show", "Enalyzer", function(object){
  options(warn = -1)
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  cat("Outlier Detection for Multiple Replicative High-throughput Data\n\n")
  cat(" Algorithm: ")

  switch(object@type, 
        Zscore = cat("Z-score criterion"),
        iqr = cat("Interquartile range (IQR) criterion"),
        siqr = cat("Semiinterquartile range (SIQR) criterion"),
        dixon = cat("Dixon's Q-test"),
        grubbs = cat("Grubbs test"),
        pair = cat("Pairwise OutlierD algoirthm"),
        diff = cat("Difference-based OutlierD algorithm"),
        proj = cat("Projection-based OutlierD algorithm")
    )
  cat(" (",object@type,")","\n")

  if( object@type %in% c("pair","diff","proj")){
    cat(" Method: ")
    switch(object@method, 
           constant = cat("constant quantile regression"),
           linear = cat("linear quantile regression"),
           nonlin = cat("non-linear quantile regression"),
           nonpar = cat("non-parametric quantile regression")
    )
    cat(" (",object@method,")","\n")
  }

  if( object@type %in% c("siqr", "iqr", "pair","diff","proj")) {
    if (object@type == "siqr") cat(" k: ", object@k, "for 2k * SIQR \n")
    else cat(" k: ",object@k,"for k * IQR\n")

    if (object@type %in% c("pair", "diff", "proj")) {
        cat(" Upper Quantile: ", object@contrl.para$Upper,"\n")
        cat(" Lower Quantile: ",object@contrl.para$Lower, "\n")
    }
  } else if(object@type == "Zscore"){
      cat(" k: ", object@k, "for |Z| > k \n")
  }

  cat(" Number of Observations: ", nrow(object@res),"\n")
  cat(" Number of Outliers: ", object@n.outliers,"\n\n")
            
  cat(" Centering: ",object@contrl.para$centering,"\n")
  cat(" Transformation: ", object@contrl.para$trans ,"\n")
            
  cat("\n Head of the Input Data \n" )
  print(object@raw.data[1:6,], digits = 3)
  cat("To see the full information of the input dataset, use a command, 'input(your_object_name)'. \n")

  cat("\n Head of the Output Results \n")

  if(object@type%in% c("Zscore", "grubbs")) out = round(object@res[, -c(1:(ncol(object@raw.data)+1))],3)
  else out = round(object@res[, -1],3) 

  out[is.na(out)] <- "."
  print(out[1:6,], digits = 3)
  cat("To see the full information for the result, use a command, 'output(your_object_name)'. \n")            

  cat('\n Head of the Peptide Numbers suspected to be an Outlier\n')
  out = rownames(object@res[object@res$Outlier,])
  print(head(out, 10))
  cat("To see the full information for the candidate outliers, use a command, 'outliers(your_object_name)'. \n")
}
)

