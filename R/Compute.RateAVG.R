#' Compute.RateAVG: function to compute the average rate parameter within an interval of the gamma distribution with parameter alpha
#'
#' This function returns a list of each step in the MCMC sampling chain
#' @param k.classes Number of discrete categories in the gamma distribution
#' @param alpha.param Alpha parameter of the gamma distribution
#' @keywords population genetics, coalescent models
#' @export
#' 
#' 

#######################################################################################################
# Function to compute the average rate values with the gamma distribution of a given set of K classes #
#######################################################################################################

Compute.RateAVG.R <- function(k.classes, alpha.param){
	quantiles.to.get <- (1:k.classes-1)/k.classes # Set up quantiles to extract
	boundaries.to.get <- rep(NA, k.classes+1) # Init vector to collect boundaries
	boundaries.to.get[1] <- 0 
	boundaries.to.get[k.classes+1] <- Inf
	mean.rates <- rep(NA, k.classes)
	
	for (i in 2:(k.classes)){
		k.quant <- quantiles.to.get[i] # Extract quantile
		k.bound <- qgamma(p = k.quant, shape = alpha.param, rate = alpha.param)
		boundaries.to.get[i] <- k.bound # Append
	}
	
	for (i in 1:k.classes){
		a = boundaries.to.get[i]
		b = boundaries.to.get[i+1]
		mean.rates[i] = pgamma(b*alpha.param, alpha.param+1, lower=T)*k.classes-pgamma(a*alpha.param, alpha.param+1, lower=T)*k.classes # compute and append mean.rates
	}
	return(mean.rates)
}
