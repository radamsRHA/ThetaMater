#' Coal.Theta.Sim: function to simulate a vector (k.vec) of mutation counts given a vector of locus lengths (l.vec), vector of sample numbers (n.vec)
#'
#' This function returns k.vec, l.vec, n.vec from simulating under the coalescent model, with the results written to 'out.file'
#' @param theta Population size parameter theta
#' @param n.list List containing the number of samples for each locus
#' @param l.list List containing the number of samples for each locus
#' @param num.loci Number of loci to simulate
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' Sim.Results <- Coal.Sim(theta = 0.01, n.vec = c(10,10,10), l.vec = c(100,100,100), num.loci = 3)
#' 
#' Returns:
#' k.vec = # of mutations for each unique data pattern
#' n.vec = # of samples for each unique data pattern
#' l.vec = # of sites for each unique data pattern
#' c.vec = # of counts of each unique data pattern


#############################################################
#  			Simulation under the coalescent mode 								#
#############################################################

Coal.Theta.Sim <- function(theta, n.vec, l.vec, num.loci){
	sim.k.vec <- rep(NA, num.loci = num.loci) # Iniatalize k vector
	sim.l.vec <- rep(NA, num.loci = num.loci) # Iniatalize k vector
	sim.n.vec <- rep(NA, num.loci = num.loci) # Iniatalize k vector
	
	n.sample <- sample(x = n.vec, size = num.loci, replace = T) # Random sample for the number of alleles per population
	l.sample <- sample(x = l.vec, size = num.loci, replace = T) # Random sample for the length of the loci
	
	for (i in 1:num.loci){ # For each simulated locus
		mutation.count <- 0 # Init mutation counter for locus
		n <- n.sample[i] # Get number of alleles to sample
		l <- l.sample[i] # Get length of locus
		k <- n # Convert to k
		
		while (k > 1) { # While 2 or more lineages to coalesce
			p.mut = (l*theta)/(k-1+(l*theta)) # Prob. of a mutation
			p.coal = (k-1)/(k-1+(l*theta)) # Prob of coalescene
			event = sample(c("M", "C"), size = 1, prob = c(p.mut, p.coal)) # Sample ranodom event
			
			if (event == "C"){ # If event is coalescence
				k = k-1 # Reduce number of lineages by 1
			}else{mutation.count = mutation.count + 1}
		}
		sim.k.vec[i] <- mutation.count # Append count of mutations
		sim.l.vec[i] <- l # Append locus length
		sim.n.vec[i] <- n # Append number of sampled alleles
	}	
	return(list(sim.k.vec = sim.k.vec, sim.n.vec = sim.n.vec, sim.l.vec = sim.l.vec)) # Return the number of mutations simulated
}
