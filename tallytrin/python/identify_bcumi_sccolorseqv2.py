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
parser.add_argument("--cmimode", default=0, type=str,
                    help='rum in CMI mode or not')

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #



log =  iotools.open_file(args.outfile  + ".log","w")


# generate set of barcodes for whitelist
barcodes = []

def find_substring(long_string, subset_str):

    pattern = re.compile(subset_str)
    match = pattern.search(long_string)

    if match:
        return match.start()
    else:
        return None


outfile = open(args.outfile, "w")

n = 0
y = 0
with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        n += 1

        seq = record.sequence
        first = 0
        for a, b in zip(Bio.pairwise2.align.localms(seq,"TCTTCCGATCT", 2, -1, -1, -1), Bio.pairwise2.align.localms(seq,"TTTTTTTTTT", 2, -1, -1, -1)):
            first +=1
            if first == 1:
                al1_a, al2_a, score_a, begin_a, end_a = a 
                al1_a, al2_b, score_b, begin_b, end_b = b 
            else:
                first = 0
                break
            length_umibarcode = len(seq[end_a:begin_b])
            barcodeumi = seq[end_a:begin_b]
            if length_umibarcode > 18 and length_umibarcode < 50:
                
                umi_start = find_substring(barcodeumi, "AGC")

                if umi_start is not None:
                    barcode = seq[end_a:end_a+10]
                    barcodes.append(barcode)
                    umi = barcodeumi[umi_start+3:umi_start+19]
                    if umi is None:
                        break
                    if len(umi) != 16:
                        break
                if args.cmimode == '1':
                    umi = seq[begin_a:end_a]
                    umi = umi[:8]
                else:
                    break
                    
                    

                seq_new = seq[begin_b:]
                quality_new = record.quality[begin_b:]
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
