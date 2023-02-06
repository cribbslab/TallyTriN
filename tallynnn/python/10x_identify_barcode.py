import sys
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
                    help='name of whitelist file')
parser.add_argument("--infile", default=None, type=str,
                    help='nanopore infile fastq  file')
parser.add_argument("--outfile", default=None, type=str,
                    help='name for output fastq files')
parser.add_argument("--cmimode", default=None, type=str,
                    help='Run the script in cmi mode for accuracy evaluation')

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #

log =  iotools.open_file(args.outfile  + ".log","w")

# generate set of barcodes for whitelist
barcodes = []

outfile = open(args.outfile, "w")

n = 0
y = 0

with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        n+=1

        first = 0

        seq = record.sequence

        for a, b in zip(Bio.pairwise2.align.localms(seq,"AGATCGGAAGAGCGT", 2, -1, -0.5, -0.1), Bio.pairwise2.align.localms(seq,"AAAAAAAAA", 2, -1, -0.5, -0.1)):
            first +=1
            if first == 1:
                al1_a, al2_a, score_a, begin_a, end_a = a 
                al1_a, al2_b, score_b, begin_b, end_b = b 
            else:
                first = 0
                break
            length_umibarcode = len(record.sequence[end_b:begin_a])

            if length_umibarcode >= 28:
                y+=1
                if args.cmimode == '1':
                    umi = record.sequence[begin_a:end_a]
                    umi = umi[:15]
                else:
                    umi = record.sequence[end_b:end_b+12]
                barcode = record.sequence[end_b+12:end_b+28]
                barcodes.append(barcode)
                seq_new = record.sequence[:begin_b]
                quality_new = record.quality[:begin_b]

                outfile.write("@%s\n%s\n+\n%s\n" % (record.name + "_" + barcode + "_" + umi, seq_new, quality_new))
                
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
