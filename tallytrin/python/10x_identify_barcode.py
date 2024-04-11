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
parser.add_argument("--barcode", default='16', type=int,
                    help='Length of barcode required')
parser.add_argument("--umi", default='12', type=int,
                    help='Length of umi required')

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


# Construct a regex pattern that matches the sequence with up to 2 mismatches in the specified middle part
sequence_pattern = r"AGATCGGAAGAGCGT"



with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        n+=1

        seq = record.sequence

        match = regex.search(sequence_pattern, seq)
        if match:
            y += 1
            barcode = seq[match.start()-16:match.start()]
            umi = seq[match.start() - 16 - args.umi:match.start() - 16]
            barcodes.append(barcode)

            seq_new = seq[:match.start()-28]
            quality_new = record.quality[:match.start()-28]

            outfile.write(f"@{record.name}_{barcode}_{umi}\n{seq_new}\n+\n{quality_new}\n")
        
       
outfile.close()

                
# Write out a list of whitelist barcodes
out_barcodes = open(args.whitelist,"w")

for i in set(barcodes):
    out_barcodes.write("%s\n" % (i))
out_barcodes.close()
        

log.close()

