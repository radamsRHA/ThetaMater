#' Read.FastaDir: Function to read a set of locus-specific fasta files provided in a directory 
#'
#' This function returns a list containing (1) the number of mutations for each RAD locus, (2) The length of each locus, and (3) the number of sequences at each locus
#' @param fasta.dir Directory containing fasta files for each individual locus
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' example.dat <- Read.FastaDir(fasta.dir = '~/Desktop/ExampleFastaDir/')

###############################################################################################################################
# Function to in a set of fasta alignments provided in a directory
###############################################################################################################################

Read.FastaDir <- function(fasta.dir){
	fasta.files <- list.files(path = fasta.dir, full.names = T)
	results.vec <- rep(NA, length(fasta.files))
	
	# loop through fasta files in directory
	for (i in 1:length(fasta.files)){
		infile.handle <- read.FASTA(file = fasta.files[i])
		locus.dat <- as.alignment(as.character(infile.handle))
		locus.seq.set <- locus.dat$seq
		n.locus <- length(locus.seq.set)
		l.locus <- length(strsplit(x = locus.seq.set[1], split = "")[[1]])
		locus.length.prime <- 0
		k.locus <- 0
		
		locus.alignment <- matrix(nrow = n.locus, ncol = l.locus)
		
		for (j in 1:n.locus){
			locus.alignment[j,] <- strsplit(x = locus.seq.set[j], split = '')[[1]]
		}
		
		for (site in 1:l.locus){
			site.aln <- locus.alignment[,site] # Extract alignment for site at locus
			if ((length(grep(pattern = 'N', x = site.aln, ignore.case = T)) == 0) && (length(grep(pattern = '-', x = site.aln)) == 0) && (length(unique(site.aln)) == 2)) { 
				k.locus <- k.locus + 1 # Add to mutation counter if biallelic site
				locus.length.prime = locus.length.prime +1
			}
			if (length(unique(site.aln))==1 && (length(grep(pattern = "N", x = site.aln, ignore.case = T)) == 0) && (length(grep(pattern = "-", x = site.aln)) == 0)){
				locus.length.prime  = locus.length.prime + 1
			}
		}
		results.vec[i] <- paste(k.locus, locus.length.prime, n.locus, sep = ",")
	}
	unique.loci <- table(results.vec)
	num.c <- length(unique.loci)
	k.vec <- rep(NA, num.c)
	l.vec <- rep(NA, num.c)
	n.vec <- rep(NA, num.c)
	c.vec <- rep(NA, num.c)
	
	for (i in 1:num.c){
		c.vec[i] <- as.numeric(unique.loci[i])
		unique.dat.0 <- names(unique.loci[i])
		unique.dat.1 <- strsplit(x = unique.dat.0, split = ",")[[1]]
		k.vec[i] <- as.numeric(unique.dat.1[1])
		l.vec[i] <- as.numeric(unique.dat.1[2])
		n.vec[i] <- as.numeric(unique.dat.1[3])
	}
	
	return(list(k.vec = k.vec, l.vec = l.vec, n.vec = n.vec, c.vec = c.vec))
}