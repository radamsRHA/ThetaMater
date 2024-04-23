#' ThetaMater.M3: (Bayesian Model 2) MCMC function that uses the MCMCmetrop1R function from MCMCpack with a discrete gamma rate variation model
#'
#' This function returns a list of each step in the MCMC sampling chain
#' @param k.vec Vector of mutation counts
#' @param l.vec Vector of locus lengths
#' @param n.vec Vector of sample numbers
#' @param c.vec Vector of data pattern counts
#' @param ngens Number of generations to run the MCMC
#' @param thin Number of generations to thin in the MCMC chain
#' @param burnin Number of generations to discard from MCMC chain
#' @param theta.shape Shape parameter of the gamma distribution for setting the prior on theta
#' @param theta.scale Scale parameter of the gamma distribution for setting the prior on theta
#' @param alpha.shape Shape parameter of the gamma distribution for setting the prior on alpha
#' @param alpha.shape Scale parameter of the gamma distribution for setting the prior on alpha
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' library(Rcpp)
#' library(ThetaMater)
#' library(MCMCpack)
#' sim.results <- Coal.Theta.Sim.File.G(theta = 0.001, n.vec = c(10,10), l.vec = c(2000,2000), num.loci = 1000, out.file = '~/Desktop/example.alleles', alpha.param = 0.5)
#' example.data <- Read.AllelesFile.NoThreshold(alleles.file = '~/Desktop/example.alleles', log.file = '~/Desktop/log.txt')
#' Theta.Posterior.G <- MCMC.Theta.M3(k.vec = example.data$k.vec, l.vec = example.data$l.vec, n.vec = example.data$n.vec, c.vec = example.data$c.vec, 
#' ngens = 10000, thin = 10, burnin = 10000, K = 4, alpha.param = 0.5)
#' 
#' 
#' 
#' 

###################################################################################################################################################
# Function to estimate the posterior probability distribution of theta given a genomic dataset (uses MCMCmetrop1R) from MCMC.pack									#
###################################################################################################################################################


ThetaMater.M3 <- function(k.vec, l.vec, n.vec, c.vec, K, ngens, burnin, thin, theta.shape, theta.scale, alpha.shape, alpha.scale){
	theta.init = rgamma(1, shape = theta.shape, scale = theta.scale)
	alpha.init = rgamma(1, shape = alpha.shape, scale = alpha.scale)
	logThetaGAMMA <- function(params, k.vec, n.vec, l.vec, c.vec, K){
		theta = abs(params[1])
		alpha.param = abs(params[2])
		rate.vec<-discrete.gamma(alpha = alpha.param, k = K)# get rate vector
		p <-  dgamma(x = alpha.param, shape = alpha.shape, scale = alpha.scale, log = T)+dgamma(x = theta, shape = theta.shape, scale = theta.scale, log = T)+ ThetaDataLikeGamma(theta = theta, kvec = k.vec, lvec = l.vec, nvec = n.vec, cvec = c.vec, rvec = rate.vec, Kclasses = K) # Proposed likelihood
		
	}
	
	MCMC.results <- MCMCmetrop1R(fun = logThetaGAMMA, theta.init = c(theta.init, alpha.init), 
															 k.vec = k.vec, n.vec = n.vec, l.vec = l.vec, 
															 c.vec = c.vec, K = K, thin=thin, mcmc=ngens, burnin=burnin, verbose = 500, logfun=T, seed = runif(n = 1, min = 0, max = 1000000), force.samp = T)
	
	
	return(MCMC.results)
}