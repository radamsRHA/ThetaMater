import msprime

length = 248956000
rho = 2e-8
number_sequences = 10
Ne = 10000
mu = 2e-8
window_size = 100000
rad_size = 1000

simulation_results = msprime.simulate(sample_size=number_sequences, Ne=Ne, length=length, recombination_rate=rho, mutation_rate=mu)


SNP_positions = []
for i in simulation_results.variants():
	genome_position = i[0]
	genome_position = int(round(genome_position, 1))
	SNP_positions.append(genome_position)


SNP_positions = set(SNP_positions)
result_aln = "0" * length
result_aln = list(result_aln)

for x in SNP_positions:
	result_aln[x-1] = "1"


result_aln = ''.join(result_aln)
th = Ne*4*mu
output_string = "Results_Rho%s_Seq%s_WindowSize%s_Radsize%s_Theta%s" % (rho, number_sequences, window_size, rad_size, th)
output = open(output_string, 'w')
for i in range(0, length, window_size):
	start = int(i)
	end = int(i) + int(rad_size)
	rad = result_aln[start:end]
	print rad
	print len(rad)
	output.write(str(rad_size) +'\t'+str(number_sequences)+'\t'+str(rad.count('1'))+'\n')



