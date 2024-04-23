#' Read.InterleavedFasta: Function to read interleaved fasta file (loci are stacked, as with Stacks output files)
#'
#' This function returns a list containing (1) the number of mutations for each RAD locus, (2) The length of each RADseq locu, and (3) the number of sequences at each RAD locus
#' @param fasta.file Fasta file (pile up across individual loci - see Stacks output file for example)
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' example.dat <- Read.InterleavedFasta(fasta.file = '~/Desktop/example.fasta')

###############################################################################################################################
# Function to read interleaved fasta files (see stacks output files) and convert into ThetaMater (infinite sites) format									#
###############################################################################################################################

Read.InterleavedFasta <- function(fasta.file){
	# Read infile data
	read.handle <- read.FASTA(file = fasta.file) # read fasta file handle
	num.lines <- length(read.handle) # Get number of line to loop through
	
	# extract locus counter from string of first locus
	results.list <- list()
	locus.counter.0 <- strsplit(x = names(as.character(read.handle[1])), split = "_")[[1]][2]
	locus.sequence.list <- list()
	sample.count = 0
	

	# Loop through lines and extract information from fasta
	for (i in 1:num.lines){
		sequence.dat <- as.character(read.handle[i]) 
		sequence.ID <- names(sequence.dat) 
		locus.ID <- strsplit(x = sequence.ID, split = "_")[[1]][2]
		
		# If same locus as locus.counter.0
		if (locus.ID == locus.counter.0){
			sample.count = sample.count + 1
			locus.length = length(sequence.dat[[1]])
			locus.sequence.list[[sample.count]] <- sequence.dat[[1]]
			
			if (i == num.lines){
				n.locus <- sample.count
				l.locus <- locus.length
				k.locus <- 0 
				locus.alignment <- matrix(nrow = n.locus, ncol = l.locus)
				locus.length.prime <- 0
				
				for (j in 1:n.locus){
					sample.seq <- locus.sequence.list[[j]]
					locus.alignment[j,] <- sample.seq
				}
				
				# Get unique sites (columns in matrix)
#				unique.col.table <- data.frame(table(apply(locus.alignment, 2, paste, collapse = ", ")))
#				num.unique.col <- length(unique.col.table[,1])
				
				# loop through unique site patterns (col)
#				for (p in 1:num.unique.col){
#					pattern.dat <- unique.col.table[p,]
#					site.pattern.0 = pattern.dat[1]
#					site.freq = as.numeric(pattern.dat[2])
#					site.pattern.1 <- gsub("\\s", "", site.pattern.0[[1]])
#					site.pattern.2 = strsplit(x = toString(site.pattern.1[[1]]), split = ',')[[1]]
					
					# If SNP site
#					if ((length(unique(site.pattern.2)) > 1) && (length(grep(pattern = "N", x = unique(site.pattern.2), ignore.case = T)) == 0) && (length(grep(pattern = "-", x = unique(site.pattern.2))) == 0)){
#						locus.length.prime = locus.length.prime + site.freq
#						k.locus = k.locus +site.freq
#					}
#					if (length(unique(site.pattern.2)) == 1){
#						locus.length.prime = locus.length.prime + site.freq
#					}
#				}
				
				for (site in 1:l.locus){ # For site in locus alignment
					site.aln <- locus.alignment[,site] # Extract alignment for site at locus
					
					if ((length(grep(pattern = 'N', x = site.aln, ignore.case = T)) == 0) && (length(grep(pattern = '-', x = site.aln)) == 0) && (length(unique(site.aln)) == 2)) { 
						k.locus <- k.locus + 1 # Add to mutation counter if biallelic site
						locus.length.prime = locus.length.prime +1
					}
					if (length(unique(site.aln))==1){
						locus.length.prime  = locus.length.prime + 1
					}
				}
				
#				
				locus.results.vec <- c(k.locus, locus.length.prime, n.locus)
				names(locus.results.vec) <- c("k.locus", "l.locus", "n.locus")
				results.list[[locus.counter.0]] <- locus.results.vec
				
			}
			
		}else{
			# Get information from current (new) sequence/locus
			sequence.dat <- as.character(read.handle[i]) 

			# Summarize previous locus information
			n.locus <- sample.count
			l.locus <- locus.length
			k.locus <- 0 
			locus.length.prime <- 0
			
			locus.alignment <- matrix(nrow = n.locus, ncol = l.locus)
			
			for (j in 1:n.locus){
				sample.seq <- locus.sequence.list[[j]]
				locus.alignment[j,] <- sample.seq
			}
			
			# Get unique sites (columns in matrix)
#			unique.col.table <- data.frame(table(apply(locus.alignment, 2, paste, collapse = ", ")))
#			num.unique.col <- length(unique.col.table[,1])
			
			# loop through unique site patterns (col)
#			for (p in 1:num.unique.col){
#				pattern.dat <- unique.col.table[p,]
#				site.pattern.0 = pattern.dat[1]
#				site.freq = as.numeric(pattern.dat[2])
#				site.pattern.1 <- gsub("\\s", "", site.pattern.0[[1]])
#				site.pattern.2 = strsplit(x = toString(site.pattern.1[[1]]), split = ',')[[1]]
				
				# If SNP site
#				if ((length(unique(site.pattern.2)) > 1) && (length(grep(pattern = "N", x = unique(site.pattern.2), ignore.case = T)) == 0) && (length(grep(pattern = "-", x = unique(site.pattern.2))) == 0)){
#					locus.length.prime = locus.length.prime + site.freq
#					k.locus = k.locus +site.freq
#				}
#				if (length(unique(site.pattern.2)) == 1){
#					locus.length.prime = locus.length.prime + site.freq
#				}
#			}
			for (site in 1:l.locus){ # For site in locus alignment
				site.aln <- locus.alignment[,site] # Extract alignment for site at locus
				
				if ((length(grep(pattern = 'N', x = site.aln, ignore.case = T)) == 0) && (length(grep(pattern = '-', x = site.aln)) == 0) && (length(unique(site.aln)) >= 2)) { # If no GAPs or N and biallelic site
					k.locus <- k.locus + 1 # Add to mutation counter if biallelic site
					locus.length.prime = locus.length.prime +1
				}
				if (length(unique(site.aln))==1){
					locus.length.prime = locus.length.prime + 1
				}
			}
			
			
			
			locus.results.vec <- c(k.locus, locus.length.prime, n.locus)
			names(locus.results.vec) <- c("k.locus", "l.locus", "n.locus")
			results.list[[locus.counter.0]] <- locus.results.vec
			
			# clear objects
			sequence.ID <- names(sequence.dat) 
			locus.ID <- strsplit(x = sequence.ID, split = "_")[[1]][2]
			sample.count = 1
			locus.sequence.list <- list()
			locus.counter.0 <- locus.ID
			locus.sequence.list[[1]] = sequence.dat[[1]]
		}
	}
	
	num.loci <- length(results.list)
	results.vec <- rep(NA, num.loci)
	
	# loop through loci
	for (i in 1:num.loci){
		locus.info.0 <- results.list[[i]]
		locus.info.1 <- toString(locus.info.0)
		results.vec[i] <- locus.info.1
	}
	
	table.results <- table(results.vec)
	num.c <- length(table.results)
	k.vec <- rep(NA, num.c)
	l.vec <- rep(NA, num.c)
	n.vec <- rep(NA, num.c)
	c.vec <- rep(NA, num.c)
	
	for (i in 1:num.c){
		c.vec[i] <- as.numeric(table.results[i])
		c.dat.0 <-  names(table.results[i])
		c.dat.1 <- gsub("\\s", "", c.dat.0)
		c.dat.2 <- strsplit(x = c.dat.1, split = ",")[[1]]
		k.vec[i] <- as.numeric(c.dat.2[1])
		l.vec[i] <- as.numeric(c.dat.2[2])
		n.vec[i] <- as.numeric(c.dat.2[3])
		
	}
	
	return(list(k.vec = k.vec, l.vec = l.vec, n.vec = n.vec, c.vec = c.vec))
}