##########################################################################################
## to use in binning pandakiller assembled samples into seqfilter parsed tag directories############################################################################################
import os
import sys
import multiprocessing as mp
from Bio import SeqIO
import string

cores=6
pool = mp.Pool(processes=cores)

out_dir = "/home/schuyler/testing_folder/Salmonella_Assembly/Franken/test"

# assemb = open(sys.argv[1], 'rU')
assemb = open("/home/schuyler/testing_folder/Salmonella_Assembly/Franken/test.fa", 'rU')


lfiles = []
for root, dirs, files in os.walk("/home/schuyler/testing_folder/Salmonella_Assembly/Franken/3_demultiplex/parse_index/trimmed_tags"):
	file_dir = root
	for file in files:
		if file.endswith("trimmed.fastq"):
			lfiles.append(file)

dict_assemb = {}

def bin(file):
	fout = open('%s/%s_assem.fasta' %(out_dir,file.split('_')[0]), 'a')
	f = open('%s/%s' %(file_dir, file) , 'rU')
	l = []
	for sequence in SeqIO.parse(f, 'fastq'):
		name = sequence.id.split(' ')[0]
		if name in dict_assemb:
			l.append(name)
	for item in l:
		SeqIO.write(dict_assemb[item], fout, 'fasta')

reads = 1
number = 0
for sequence in SeqIO.parse(assemb, "fasta"):
	sequence.id = sequence.id.rsplit(':', 1)[0]
	sequence.description = sequence.description.rsplit(':', 1)[0]
	sequence.name = sequence.name.rsplit(':', 1)[0]
	dict_assemb[sequence.id] = sequence
	number += 1
	if number == 1000000: #200,000 reads per 1GB of
		pool.map(bin, lfiles)

		print "%s Million Reads Processed" %reads
		reads += 1
		number = 0
		dict_assemb = {}


		




