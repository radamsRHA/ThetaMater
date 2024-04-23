#' ThetaMater.M2: (Bayesian Model 2) MCMC function that uses the MCMCmetrop1R function from MCMCpack with a discrete gamma rate variation model
#'
#' This function returns a list of each step in the MCMC sampling chain
#' @param k.vec Vector of mutation counts
#' @param l.vec Vector of locus lengths
#' @param n.vec Vector of sample numbers
#' @param c.vec Vector of data pattern counts
#' @param alpha.param Shape of the gamma distribution for describing the amount of among-locus rate variation
#' @param ngens Number of generations to run the MCMC
#' @param K Number of classese to approximate the gamma distribution
#' @param thin Number of generations to thin in the MCMC chain
#' @param burnin Number of generations to discard from MCMC chain
#' @param theta.shape Shape parameter of the gamma distribution for setting the prior on theta
#' @param theta.scale Scale parameter of the gamma distribution for setting the prior on theta
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' library(Rcpp)
#' library(ThetaMater)
#' library(MCMCpack)
#' 

###################################################################################################################################################
# Function to estimate the posterior probability distribution of theta given a genomic dataset (uses MCMCmetrop1R) from MCMC.pack	and a fixed alpha value of among-locus rate variation								#
###################################################################################################################################################


ThetaMater.M2 <- function(k.vec, l.vec, n.vec, c.vec, K, alpha.param, ngens, burnin, thin, theta.shape, theta.scale){
	rate.vec<-discrete.gamma(alpha = alpha.param, k = K)# get rate vector
	#theta.init = runif(n = 1, min = 0, max = 1) # randomly sample starting position
	theta.init = rgamma(n = 1, shape = theta.shape, scale = theta.scale)
	logThetaGAMMA <- function(theta, k.vec, n.vec, l.vec, c.vec, r.vec){
		theta=abs(theta)
		p <- dgamma(x = theta, shape = theta.shape, scale = theta.scale, log = T)+ ThetaDataLikeGamma(theta = theta, kvec = k.vec, lvec = l.vec, nvec = n.vec, cvec = c.vec, rvec = r.vec, Kclasses = K) # Proposed likelihood
	}
	
	MCMC.results <- MCMCmetrop1R(fun = logThetaGAMMA, theta.init = theta.init, 
																					k.vec = k.vec, n.vec = n.vec, l.vec = l.vec, 
																					c.vec = c.vec, r.vec = rate.vec, thin=thin, mcmc=ngens, burnin=burnin, verbose = 500, logfun=T, seed = runif(n = 1, min = 0, max = 1000000))
	

	return(MCMC.results)
}