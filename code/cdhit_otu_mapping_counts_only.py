import sys
import os
import re
import numpy
from itertools import groupby
def clstr_iter(cdhit_clstr):
        f = open(cdhit_clstr)
        citer = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
        for cluster in citer:
                seq = []
                for line in citer.next():
                        if "*" in line:
                                string = re.split("\t| |>|\||;", line)
                                cluster = string[3]
                        else: 
                                string = re.split("\t| |;|/|>|\.|\||,", line)
                                read_count = int(string[5].split("=")[1])
                                seq.append(read_count)
                yield cluster, seq
clstr = sys.argv[1]
d = dict(clstr_iter(clstr))
for item in d:
        print "%s\t%s\t%s" % (item.strip("."), clstr.rsplit(".")[0].rsplit("/")[-1].rsplit("_")[-1], sum(d[item]))
