#' FilterData.PPS: function to filter data using a threshold provided by ThetaMater.PPS
#'
#' This function filters a given dataset using a threshold of mutaiton counts. All loci with mutation acounts above this threshold will be removed from analysis
#' @param threshold Maximum number of mutations used to filter loci
#' @param dataset Target dataset to be filtered (in infinite sites format)
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' library(Rcpp)
#' library(ThetaMater)
#' library(MCMCpack)


###################################################################################################################################################
# Function to filter datasets given a threshold inferred via PPS
###################################################################################################################################################


FilterData.PPS <- function(dataset, threshold){
	num.unique.loci <- length(dataset$c.vec)
	counter = 0
	for (i in 1:num.unique.loci){
		k.locus <- dataset$k.vec[i]
		if (k.locus<= threshold){
			counter = counter + 1
		}
	}
	
	k.vec <- rep(NA, counter)
	l.vec <- rep(NA, counter)
	n.vec <- rep(NA, counter)
	c.vec <- rep(NA, counter)
	counter.2 = 0 
	
	for (j in 1:num.unique.loci){
		k.locus <- dataset$k.vec[j]
		l.locus <- dataset$l.vec[j]
		n.locus <- dataset$n.vec[j]
		c.locus <- dataset$c.vec[j]
		if (k.locus <= threshold){
			counter.2 = counter.2 + 1
			k.vec[counter.2] <- k.locus
			l.vec[counter.2] <- l.locus
			n.vec[counter.2] <- n.locus
			c.vec[counter.2] <- c.locus
		}
	}
	return(list(k.vec = k.vec, l.vec = l.vec, n.vec = n.vec, c.vec = c.vec))
}