## This script takes Rachel's fasta reads with headers
# >barcodelabel=B1000;_M02357:0:000000000-AA7V9:1:1101:17867:3631;size=259435;
# and splits into fasta files by barcode= " " and makes new headers as
# >M02357:0:000000000-AA7V9:1:1101:17867 size=259435


import os
import sys
from Bio import SeqIO
import string


out_dir = sys.argv[2]
assemb = open(sys.argv[1], 'rU')

for sequence in SeqIO.parse(assemb, "fasta"):
        barcode = sequence.description.rsplit(':', 1)[0].split(";")[0].split("=")[1]
        sequence.description = sequence.description.rsplit(';', 1)[0].split(";")[-1]
        sequence.id = sequence.id.rsplit(':', 1)[0].split("_")[1] + ";" + sequence.description
        sequence.description = ""
        sequence.name = sequence.name.rsplit(':', 1)[0].split("_")[1]
        fout = open('%s/%s.fasta' %(out_dir,barcode.split('.fasta')[0]), 'a')
        SeqIO.write(sequence, fout, 'fasta')
