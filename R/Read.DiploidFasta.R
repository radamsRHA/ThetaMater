#' Read.DiploidFasta: Function to read diploid genome sequence in fasta format
#'
#' This function returns a list containing (1) the number of mutations for each RAD locus, (2) The length of each RADseq locu, and (3) the number of sequences at each RAD locus
#' @param genome.fasta.file Fasta file comprising multiple contigs with ambiquity codes 
#' @keywords population genetics, coalescent models
#' @export
#' @examples
#' Read.DiploidFasta(genome.fasta.file = '~/Desktop/example.genome.fasta', num.contigs = 4)

###############################################################################################################################
# Function to read diploid (SNPs coded as ambiquities) genome fasta file
###############################################################################################################################

Read.DiploidFasta <- function(genome.fasta.file, num.contigs){
	infile.handle <- read.FASTA(file = genome.fasta.file)
	total.contigs <- length(infile.handle)
	
	# Randomly sample contigs for thetamater
	sample.random.contigs <- sample(x = names(infile.handle), size = num.contigs, replace = F)
	results.vec <- rep(NA, num.contigs)
	
	for (i in 1:num.contigs){
		contig.header <- sample.random.contigs[i]
		contig.data <- as.alignment(as.character(infile.handle[contig.header]))
		contig.seq <- contig.data$seq
		count.N <- str_count(string = contig.seq, pattern = ignore.case("N"))
		count.gap <- length(grep(pattern = '-',x = contig.seq, ignore.case = T))
		contig.items <- nchar(contig.seq)
		contig.length = contig.items-(count.N+count.gap)
		
		A.count <- str_count(string = contig.seq, pattern = ignore.case("A"))
		C.count <- str_count(string = contig.seq, pattern = ignore.case("C"))
		G.count <- str_count(string = contig.seq, pattern = ignore.case("G"))
		T.count <- str_count(string = contig.seq, pattern = ignore.case("T"))
		k <- contig.length-(A.count+C.count+G.count+T.count)
		results.vec[i] <- paste(k, contig.length, sep = " ")
		
	}
	
	unique.contig.dat <- table(results.vec)
	num.unique.contigs <- length(unique.contig.dat)
	k.vec <- rep(NA, num.unique.contigs)
	n.vec <- rep(2, num.unique.contigs)
	l.vec <- rep(NA, num.unique.contigs)
	c.vec <- rep(NA, num.unique.contigs)
	
	
	for (i in 1:num.unique.contigs){
		c.vec[i] <- as.numeric(unique.contig.dat[i])
		
		unique.dat.split <- strsplit(names(unique.contig.dat[i]), split = " ")[[1]]
		k.vec[i] <- as.numeric(unique.dat.split[1])
		l.vec[i] <- as.numeric(unique.dat.split[2])
		
	}
	return(list(k.vec = k.vec, l.vec = l.vec, n.vec = n.vec, c.vec = c.vec))
}