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

#gc.disable()
n = 0
for bed1, bed2 in zip(bed1, bed2):

    barcode = bed1.split("\t")[4]
    ig = bed1.split("\t")[3]
    transloc = bed2.split("\t")[3]
    trans_gene = ig + "_" + transloc + "_"+ barcode.strip()

    total += 1

    if  "ENSG" in trans_gene and "ENSM" in trans_gene:
        mouse_human +=1
    trans_list.append(str(trans_gene))

out_table.write("gene1\tgene2\ttrans\tbarcode\tcount\n")
out_log.write("Total = %s\n mouse and human fusions = %s"%(total, mouse_human))

trans_counter =collections.Counter(trans_list)

for k in trans_counter:
    ig, gene, barcode = k.split("_")
    if ig == gene:
        pass
    else:
        gene2 = ig + "_" + gene
        counter = trans_counter[k]

        out_table.write("%s\t%s\t%s\t%s\t%s\n"%(ig, gene, gene2, barcode, str(counter)))
    
out_table.close()
