import sys
import re
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import Bio
import Bio.pairwise2
from Bio.pairwise2 import format_alignment

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("identify_perfect_nano.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--whitelist", default=None, type=str,
                    help='a file naming the outfile for the whitelist of barcodes.')
parser.add_argument("--infile", default=None, type=str,
                    help='nanopore infile fastq  file')
parser.add_argument("--outfile", default=None, type=str,
                    help='name for output fastq files')
parser.add_argument("--smarter", default=None, type=str,
                    help='specify if smarter UMI sequence is used')
parser.add_argument("--trimer", default=None, type=str,
                    help='specify if the trimer UMI or normal UMI should be used')



args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #



log =  iotools.open_file(args.outfile  + ".log","w")


# generate set of barcodes for whitelist
barcodes = []

def poly_a_length(input_string):
    ''' find ployA tail and then calculate length and return the end A string index)'''

    poly_a_pattern = re.compile(r"A+")

    poly_a_matches = re.findall(poly_a_pattern, input_string)

    longest_poly_a = max(poly_a_matches, key=len)
    longest_poly_a_length = len(longest_poly_a)
    last_a_index = input_string.index(longest_poly_a) + longest_poly_a_length

    return last_a_index


def most_common(lst):
    return max(set(lst), key=lst.count)

outfile = open(args.outfile, "w")

n = 0
y = 0


if smarter:

    with pysam.FastxFile(args.infile) as fh:
    
        for record in fh:
            n += 1

            seq = record.sequence
            first = 0
            for a, b in zip(Bio.pairwise2.align.localms(seq,"GTACTCTGCGTT", 2, -1, -1, -1), Bio.pairwise2.align.localms(seq,"AAAAAAAAA", 2, -1, -1, -1)):
                first +=1
                if first == 1:
                    al1_a, al2_a, score_a, begin_a, end_a = a 
                    al1_a, al2_b, score_b, begin_b, end_b = b 
                else:
                    first = 0
                    break
                    length_umibarcode = len(seq[end_b:begin_a])

                if length_umibarcode > 48:
                
                    barcode = seq[begin_a -30:begin_a]
                    barcodes.append(barcode)
                    umi = seq[end_b:end_b + 21]
                    seq_new = seq[:begin_b]
                    quality_new = record.quality[:begin_b]
                    y += 1
                    outfile.write("@%s\n%s\n+\n%s\n" % (record.name + "_" + barcode + "_" + umi, seq_new, quality_new))
                else:
                    pass
                
    outfile.close()

else:
    with pysam.FastxFile(args.infile) as fh:
    
        for record in fh:
            n += 1

            seq = record.sequence
            first = 0
            for a, b in zip(Bio.pairwise2.align.localms(seq,"GTACTCTGCGTT", 2, -1, -1, -1), Bio.pairwise2.align.localms(seq,"AAAAAAAAA", 2, -1, -1, -1)):
                first +=1
                if first == 1:
                    al1_a, al2_a, score_a, begin_a, end_a = a 
                    al1_a, al2_b, score_b, begin_b, end_b = b 
                else:
                    first = 0
                    break
                    length_umibarcode = len(seq[end_b:begin_a])

                if length_umibarcode > 48:
                
                    barcode = seq[begin_a -30:begin_a]
                    barcodes.append(barcode)
                    umi = seq[end_b:end_b + 21]
                    seq_new = seq[:begin_b]
                    quality_new = record.quality[:begin_b]
                    y += 1
                    outfile.write("@%s\n%s\n+\n%s\n" % (record.name + "_" + barcode + "_" + umi, seq_new, quality_new))
                else:
                    pass
                
    outfile.close()

                
# Write out a list of whitelist barcodes
out_barcodes = open(args.whitelist,"w")

for i in set(barcodes):
    out_barcodes.write("%s\n" % (i))
out_barcodes.close()
        

log.write("The number of total reads: %s\n" %(n))
log.write("The number of total reads that have a correct barcode and UMI: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((y/n)*100))

log.close()
