import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import re
import collections



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("complement_ployA.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--outfile", default=None, type=str,
                    help='ouput counts  file')
parser.add_argument("--bed1", default=None, type=str,
                    help='bedfile1 out  file')
parser.add_argument("--bed2", default=None, type=str,
                    help='bedfile2 out  file')

args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #

bed1 = open(args.bed1, "r")
bed2 = open(args.bed2, "r")

out_table = open(args.outfile, "w")
out_log = open(args.outfile + ".log", "w")

trans_list = []


mouse_human = 0
total = 0

trans_list = []

mouse_human = 0
total = 0
n = 0
for bed1, bed2 in zip(bed1, bed2):


    umi = bed1.split("\t")[4]
    ig = bed1.split("\t")[3]
    transloc = bed2.split("\t")[3]
    trans_gene = ig + "_" + transloc 

    total += 1

    if  "ENSG" in trans_gene and "ENSM" in trans_gene:
        mouse_human +=1
    trans_list.append(str(trans_gene))

out_table.write("trans\tcount\n")
out_log.write("Total = %s\n mouse and human fusions = %s"%(total, mouse_human))


trans_counter =collections.Counter(trans_list)

for k in trans_counter:
    ig, gene = k.split("_")
    if ig == gene:
        pass
    else:
        gene2 = ig + "_" + gene
        counter = trans_counter[k]
        print(gene2, counter)

        out_table.write("%s\t%s\n"%(gene2, str(counter)))
    
out_table.close()
