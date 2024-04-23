#' Coal.Theta.Sim.File: function to simulate a vector (k.vec) of mutation counts given a vector of locus lengths (l.vec), vector of sample numbers (n.vec)
#'
#' This function returns k.vec, l.vec, n.vec from simulating under the coalescent model, with the results written to 'out.file'
#' @param theta Population size parameter theta
#' @param n.list List containing the number of samples for each locus
#' @param l.list List containing the number of samples for each locus
#' @param num.loci Number of loci to simulate
#' @param out.file Outfile to write in alleles format
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' 
#' example.data <- Coal.Theta.Sim.File(theta = 0.017, n.vec = c(5,7), l.vec = c(1,1), num.loci = 100, out.file = '~/Desktop/example.alleles')
#' 
#' Coal.Theta.Sim.File Returns:
#' k.vec = # of mutations for each unique data pattern
#' n.vec = # of samples for each unique data pattern
#' l.vec = # of sites for each unique data pattern
#' c.vec = # of counts of each unique data pattern


#############################################################
#  			Simulation under the coalescent mode 								#
#############################################################

Coal.Theta.Sim.File <- function(theta, n.vec, l.vec, num.loci, out.file){
	nucleotide.set <- c("A", "C", "G", "T") # Set of nucleotides
	if (file.exists(out.file)) file.remove(out.file) # Clear output
	
	sim.k.vec <- rep(NA, num.loci = num.loci) # Initialize k vector
	sim.l.vec <- rep(NA, num.loci = num.loci) # Initialize k vector
	sim.n.vec <- rep(NA, num.loci = num.loci) # Initialize k vector
	i <- 1 # Initialize increment counter
	locus.counter <- 0 # Initialize 'good' locus counter, number of mutations less than length of locus (infinite sites model)
	
	n.sample <- sample(x = n.vec, size = num.loci, replace = T) # Random sample for the number of alleles per population
	l.sample <- sample(x = l.vec, size = num.loci, replace = T) # Random sample for the length of the loci

	while (i <= num.loci){ # For each simulated locus
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
		
		if (mutation.count <= l){
			i = i + 1 # increase locus counter
			locus.counter = locus.counter + 1 # Increase 'good' locus counter (k < locus.length; infinite sites model)
			sim.k.vec[locus.counter] <- mutation.count # Append count of mutations
			sim.l.vec[locus.counter] <- l # Append locus length
			sim.n.vec[locus.counter] <- n # Append number of sampled alleles
			null.string <- sample(x = nucleotide.set, size = l, replace = T) # Sample ancestral sequence without mutations
		
			for (j in 1:n){ # For allele sample
				string1 <- gsub(pattern = 'XXX', replacement = j, x = ">SimLocus_YYY_Sample_XXX  ") # Replace
				string2 <- gsub(pattern = 'YYY', replacement = locus.counter, x = string1) # Replace
				sequence.string <- null.string # Assign ancestral sequence
			
				if (j == 1 && mutation.count > 0){ # If at the first sample
					for (sim.k in 1:mutation.count){ # For each simulated mutation
						ancestral.allele <- sequence.string[sim.k] # Get ancestral allele
						derived.nuc.set <- nucleotide.set[nucleotide.set != ancestral.allele] # Reduced set for mutation
						derived.allele <- sample(x = derived.nuc.set, size = 1, replace = T) # Sample mutated allele
						sequence.string[sim.k] <- derived.allele # Assign derived allele
					}
				}
				string3 <- gsub(pattern = 'ZZZ', replacement = string2, x = "ZZZ  WWW")
				string4 <- gsub(pattern = 'WWW', replacement = toString(paste(sequence.string, sep = '', collapse = '')), x = string3)
				write(x = string4, file = out.file, append = T) # Write results to outfile
			}
			print(locus.counter)
			write(x = "//                                                                                             |1", file = out.file, append = T)
		}
	}
	return(list(sim.k.vec = sim.k.vec, sim.n.vec = sim.n.vec, sim.l.vec = sim.l.vec)) # Return the number of mutations simulated
}
