import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
from collections import Counter


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("correct_10xbarcode.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--outfile", default=None, type=str,
                    help='name for output fastq files')
parser.add_argument("--whitelist", default=None, type=str,
                    help='The combined whitelist barcode file')
parser.add_argument("--cells", default=None, type=str,
                    help='The combined whitelist barcode file')
parser.add_argument("--cmimode", default=None, type=str,
                    help='Specify if CMI mode is activated')
parser.add_argument("--umi", default=None, type=str,
                    help='UMI lengths')
args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #


outf =  iotools.open_file(args.outfile,"w")
log =  iotools.open_file(args.outfile + ".log","w")
barcodes = iotools.open_file(args.whitelist)

barcodes_new = []

for line in barcodes:
    bc = line.strip()
    barcodes_new.append(bc)

lists = Counter(barcodes_new).most_common(int(args.cells))
barcode_lists = []

for key, count in lists:
    barcode_lists.append(key)


def closest_match(input_string, string_list):
    closest_match = None
    min_distance = float('inf')

    for candidate in string_list:
        distance = sum(ch1 != ch2 for ch1, ch2 in zip(input_string, candidate))
        if distance <= 2 and distance < min_distance:
            closest_match = candidate
            min_distance = distance
    return closest_match


with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        
        seq = record.sequence
        name = record.name.split("_")[0]
        barcode = record.name.split("_")[1]
        umi = record.name.split("_")[2]

        match = closest_match(barcode, barcodes)
        
        if match:
            barcode = match
        else:
            pass    
        if len(umi) == args.umi:
            outf.write("@%s\n%s\n+\n%s\n" % (name + "_" + barcode + "_" + umi, record.sequence, record.quality))
            

outf.close()
barcodes.close()
