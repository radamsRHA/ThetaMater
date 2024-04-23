#' ThetaMater.M1: (Bayesian Model 2) MCMC function that uses the MCMCmetrop1R function from MCMCpack
#'
#' This function returns a list of each step in the MCMC sampling chain
#' @param k.vec Vector of mutation counts
#' @param l.vec Vector of locus lengths
#' @param n.vec Vector of sample numbers
#' @param c.vec Vector of data pattern counts
#' @param ngens Number of generations to run the MCMC
#' @param burnin Number of generations to discard as burnin
#' @param thin Number of generations to thin in the MCMC chain
#' @param theta.shape Shape parameter of the gamma distribution for setting the prior on theta
#' @param theta.scale Scale parameter of the gamma distribution for setting the prior on theta
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' library(Rcpp)
#' library(ThetaMater)
#' library(MCMCpack)


###################################################################################################################################################
# Function to estimate the posterior probability distribution of theta given a genomic dataset (uses MCMCmetrop1R) from MCMC.pack									#
###################################################################################################################################################


ThetaMater.M1 <- function(k.vec, l.vec, n.vec, c.vec, ngens, burnin, thin, theta.shape, theta.scale){
	#theta.init = runif(n = 1, min = theta.min, max = theta.max) # init sample of theta
	theta.init = rgamma(n = 1, shape = theta.shape, scale = theta.scale)
	logThetaConstant <- function(theta, k.vec, n.vec, l.vec, c.vec){
		theta = abs(theta)
		#p <- dunif(theta, min = theta.min, max = theta.max, log = T)+ ThetaDataLike(theta = theta, kvec = k.vec, lvec = l.vec, nvec = n.vec, cvec = c.vec) # Proposed likelihood
		p <- dgamma(x = theta, shape = theta.shape, scale = theta.scale, log = T)+ ThetaDataLike(theta = theta, kvec = k.vec, lvec = l.vec, nvec = n.vec, cvec = c.vec) # Proposed likelihood
	}
		MCMC.results <- MCMCmetrop1R(fun = logThetaConstant, theta.init = c(theta.init), 
																							 k.vec = k.vec, n.vec = n.vec, l.vec = l.vec, 
																							 c.vec = c.vec,  thin=thin, mcmc=ngens, burnin=burnin, verbose = 500, logfun=T, seed = runif(n = 1, min = 0, max = 1000000000))
		
	return(MCMC.results)
}