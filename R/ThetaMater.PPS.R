#' ThetaMater.PPS: function to simulate posterior predictive distributions of mutation counts using the posterior distributions estimated using ThetaMater
#'
#' This function returns a vector of mutation counts for every simulated locus (the posterior predictive distribution)
#' @param theta.MCMC MCMC chain sampling of theta parameter provided by ThetaMater
#' @param l.vec Vector of locus lengths
#' @param n.vec Vector of sample numbers
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' sim.results <- Coal.Theta.Sim.File(theta = .9/10, n.vec = c(5,5), l.vec = c(100,100), num.loci = 1000, out.file = '~/Desktop/example.alleles')
#' example.data <- Read.AllelesFile.Threshold(alleles.file = '~/Desktop/example.alleles', threshold = 1000, log.file = '~/Desktop/log.txt')
#' Theta.Posterior <- ThetaMater.M3(k.vec = example.data$k.vec, l.vec = example.data$l.vec, n.vec = example.data$n.vec, c.vec = example.data$c.vec, ngens = 10000)
#' 

#################################################################################################################
# Function to simulate a posterior distribution of mutation counts given the results of ThetaMater							#
#################################################################################################################


ThetaMater.PPS <- function(theta.MCMC, l.vec, n.vec){
	num.MCMC.steps <- length(theta.MCMC) # get number of MCMC steps in posterior distribution
	PPS.results <- rep(NA, num.MCMC.steps) # init results
	
	for (a in 1:num.MCMC.steps){
		theta<- theta.MCMC[a] # extract individual theta estimate for MCMC step
		n <- sample(x = n.vec, size = 1)
		l <- sample(x = l.vec, size = 1)
		mutation.count <- 0
		k <- n # Convert to k
		
		while (k > 1) { # While 2 or more lineages to coalesce
			p.mut = (l*theta)/(k-1+(l*theta)) # Prob. of a mutation
			p.coal = (k-1)/(k-1+(l*theta)) # Prob of coalescene
			event = sample(c("M", "C"), size = 1, prob = c(p.mut, p.coal)) # Sample ranodom event
			if (event == "C"){ # If event is coalescence
				k = k-1 # Reduce number of lineages by 1
			}else{mutation.count = mutation.count + 1}
		}
		
		
		PPS.results[a] <- mutation.count
		
		
		#PPS.results[a] <- Coal.Theta.Sim.File(theta = theta.sample, l.vec = l, n.vec = n, num.loci = 1,out.file = "temp.txt")$sim.k.vec
		
	}
	return(PPS.results)
}
