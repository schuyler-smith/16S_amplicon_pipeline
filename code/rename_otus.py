import sys
import os

with open(sys.argv[1]) as file:
	OTU = 0
	for line in file:
		if ">" in line:
			line = ">OTU_%s" %OTU
			OTU += 1
		print line.strip()
