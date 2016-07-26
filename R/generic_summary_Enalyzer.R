setGeneric("summary")
setMethod("summary", "Enalyzer", function(object){

})

summary.Enalyzer <- function(object, ...){
	cat("\nCall:\n")
	print(object@call)
	
	cat("\n Analysis of Energy Informatics Data")

}