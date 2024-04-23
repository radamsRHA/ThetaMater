#' ThetaMater.M0: (Bayesian Model 1) Simple MCMC algorithm to estimate the posteriori distribution of theta from genomic data
#'
#' This function returns a list of each step in the MCMC sampling chain
#' @param k.vec Vector of mutation counts
#' @param l.vec Vector of locus lengths
#' @param n.vec Vector of sample numbers
#' @param c.vec Vector of data pattern counts
#' @param ngens Number of generations to run the MCMC
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' sim.results <- Coal.Theta.Sim.File(theta = .9/10, n.vec = c(5,5), l.vec = c(100,100), num.loci = 1000, out.file = '~/Desktop/example.alleles')
#' example.data <- Read.AllelesFile.Threshold(alleles.file = '~/Desktop/example.alleles', threshold = 1000, log.file = '~/Desktop/log.txt')
#' Theta.Posterior <- MCMC.Theta.M1.R(k.vec = example.data$k.vec, l.vec = example.data$l.vec, n.vec = example.data$n.vec, c.vec = example.data$c.vec, ngens = 10000)
#' 

#################################################################################################################
# Function to estimate the posterior probability distribution of theta given a genomic dataset									#
#################################################################################################################


ThetaMater.M0 <- function(k.vec, l.vec, n.vec, c.vec, ngens){
	chain = array(dim = c(ngens+1,1)) # Init chain to collect MCMC samples
	chain[1,] = runif(n = 1, min = 0.001, max = 0.01) # Starting value of theta is uniformly sampled
	for (i in 1:ngens){ # For each chain to be sample
		print(i)
		ProposedTheta <- abs(rnorm(1,mean = chain[i,], sd= 0.0001)) # Absolute value for the proposed theta
		ProposedProb <- dunif(ProposedTheta, min = 0.001, max = 0.01, log = T)+ ThetaDataLike(theta = ProposedTheta, kvec = k.vec, lvec = l.vec, nvec = n.vec, cvec = c.vec) # Proposed likelihood
		
		CurrentProb <- dunif(chain[i,], min = 0.001, max = 0.01, log = T) + ThetaDataLike(theta = chain[i,], kvec = k.vec, lvec = l.vec, nvec = n.vec, cvec = c.vec) # Current Likelihood
		probab <- exp(ProposedProb - CurrentProb) # Difference between the two probabilities
		if (runif(n = 1, min = 0, max = 1)<probab){ # If less than ratio of prior*likelihood odds
			chain[i+1,] = ProposedTheta # Add to proposed theta to chain
		}else{
			chain[i+1,] = chain[i,] # Add current to chain
		}
	}
	return(chain) # return MCMC chain
}
