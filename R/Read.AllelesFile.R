#' Read.AllelesFile: Function to read PyRAD *.alleles files and convert to ThetaMater format
#'
#' This function returns a list containing (1) the number of mutations for each RAD locus, (2) The length of each RADseq locu, and (3) the number of sequences at each RAD locus
#' @param allele.file pyRAD output file (*.alleles)
#' @param log.file Log file to write filtering results
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' example.dat <- Read.AllelesFile(alleles.file = '~/Desktop/example.alleles')

#######################################################################################################
# Function to read *alleles file and convert into ThetaMater (infinite sites) format									#
#######################################################################################################



Read.AllelesFile <- function(alleles.file){

	# Read alleles.file data
	read.handle <- readLines(alleles.file) # Handle to readlines from alleles file
	num.loci <- length(grep(pattern = '//', x = read.handle)) # Get total number of loci
	num.lines <- length(read.handle) # Get number of line to loop through
	
	# Initiate counters for reading alleles.file
	locus.count = 1 # Iniate with with the first locus
	sample.count = 0 # Initiate sample counter per locus (restarts on line 32)
	
	# Initiate vectors for summarize alleles.file for ThetaMater
	k1.vec <- rep(NA, num.loci) # Collect k for passing loci
	l1.vec <- rep(NA, num.loci) # Collect n for passing loci
	n1.vec <- rep(NA, num.loci) # Collect l for passing loci
	nl1.vec <- rep(NA, num.loci) # Collect combined n and l for passing loci
	

	
	print(gsub(pattern = 'YYY', replacement = num.loci, gsub(pattern = 'XXX', replacement = alleles.file, x = 'Reading XXX for all YYY loci...')))
	
	for (i in 1:num.lines) { # For each line in alleles.file
		line.string <- read.handle[i] # Extract line string
		
		if (length(grep(pattern = '//', x = line.string)) == 0) { # If not a new locus		
			sample.count = sample.count + 1 # sample.counter +=1
			seq.seq <- strsplit(x = strsplit(x= line.string, split = '\\s+')[[1]][2], split = '')[[1]] # Split line into sequence
			locus.length <- length(seq.seq) # Get length of locus
			
		}else{
			n.locus <- sample.count # Get sample count for locus
			l.locus <- locus.length # Get number of sites at locus
			
			k.locus <- 0 # Initiate locus mutation counter 
			locus.alignment <- matrix(nrow = n.locus, ncol = l.locus) # Initiate Matrix for locus alignment
			ID.vec <- rep(NA, n.locus) # Initiate vector for sequence IDs
			
			for (j in n.locus:1) { # For each sample at locus
				sample.string <- read.handle[i-j] # Get string for sample
				seq.ID <- strsplit(x= sample.string, split = '\\s+')[[1]][1] # Split line to get sequence ID
				seq.seq <- strsplit(x = strsplit(x= sample.string, split = '\\s+')[[1]][2], split = '')[[1]] # Split line into sequence
				locus.alignment[j,] <- seq.seq # Append to alignment
				ID.vec[j] <- seq.ID # Append SeqID
				
			}
			
			rownames(locus.alignment) <- ID.vec # Append SeqID names
			
			for (site in 1:l.locus){ # For site in locus alignment
				site.aln <- locus.alignment[,site] # Extract alignment for site at locus
				
				if ((length(grep(pattern = 'N', x = site.aln, ignore.case = T)) == 0) && (length(grep(pattern = '-', x = site.aln)) == 0) && (length(unique(site.aln)) == 2)) { # If no GAPs or N and biallelic site
					k.locus <- k.locus + 1 # Add to mutation counter if biallelic site
				}
			}
			
			k1.vec[locus.count] <- k.locus # Add to k 
			n1.vec[locus.count] <- n.locus # Add to n
			l1.vec[locus.count] <- l.locus # Add to l
			nl1.vec[locus.count] <- paste(n.locus, l.locus, k.locus, sep = ',') # Add to nl vector
			
			
			
			# RESTART COUNTERS #
			sample.count = 0
			locus.count = locus.count + 1 # Add to locus counter
		}
	}
	
	num.unique.patterns = length(table(nl1.vec)) # Number of unique patterns in the dataset
	k.vec <- rep(NA, num.unique.patterns) # Collect unique k count
	l.vec <- rep(NA, num.unique.patterns) # Collect unique l count
	n.vec <- rep(NA, num.unique.patterns) # Collect unique n count
	c.vec <- rep(NA, num.unique.patterns) # Collect unique pattern counts
	
	for (i in 1:num.unique.patterns){ # For each unique pattern in the dataset
		unique.pattern <- table(nl1.vec)[i] # Extract data and count
		pattern.n <- strsplit(x = names(unique.pattern), split = ',')[[1]][1] # Extract unique k
		pattern.l <- strsplit(x = names(unique.pattern), split = ',')[[1]][2] # Extract unique l
		pattern.k <- strsplit(x = names(unique.pattern), split = ',')[[1]][3] # Extract unique n
		
		k.vec[i] <- as.numeric(pattern.k) # Append k
		l.vec[i] <- as.numeric(pattern.l) # Append l
		n.vec[i] <- as.numeric(pattern.n) # Append n
		c.vec[i] <- as.numeric(unique.pattern) # Append c
	}
	
	print(gsub(pattern = 'ZZZ', replacement = num.loci, x = 'ALL DONE! returning a list of k.vec, n.vec, l.vec, locus.nums for ZZZ loci'))
	return(list(k.vec = k.vec, n.vec = n.vec, l.vec = l.vec, c.vec = c.vec))
}