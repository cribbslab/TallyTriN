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
                    help='A file naming the outfile for the whitelist of barcodes.')
parser.add_argument("--infile", default=None, type=str,
                    help='Nanopore infile FASTQ file')
parser.add_argument("--outfile", default=None, type=str,
                    help='Name for output FASTQ files')
parser.add_argument("--barcode_len", default=None, type=str,
                    help='Length of each cell barcode')
parser.add_argument("--cmimode", default=None, type=str,
                    help='Run the script in cmi mode for accuracy evaluation')

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code ############################################# #
# ########################################################################### #

log = iotools.open_file(args.outfile + ".log", "w")

# Generate set of barcodes for whitelist
barcodes = []


def extract_barcodes_and_umi(seq):
    """
    Extract barcodes 1, 2, 3, and the UMI sequence from the sequence based on the new structure.
    If any component is missing, return None.
    """
    pattern = r"AAGCAGTGGTATCAACGCAGAGTAC(?P<barcode1>.{10})CAGCTACTGC(?P<barcode2>.{10})CGAGTAGTACCCT(?P<barcode3>.{10})GATAGAGCCG(?P<umi>.{8})"
    match = re.search(pattern, seq)
    if match:
        barcode1 = match.group("barcode1")
        barcode2 = match.group("barcode2")
        barcode3 = match.group("barcode3")
        umi = match.group("umi")
        return barcode1, barcode2, barcode3, umi
    return None  # Return None if the read does not match the expected structure


outfile = open(args.outfile, "w")

n = 0
y = 0
with pysam.FastxFile(args.infile) as fh:
    for record in fh:
        n += 1
        seq = record.sequence

        # Extract barcodes and UMI
        result = extract_barcodes_and_umi(seq)

        if result is not None:
            barcode1, barcode2, barcode3, umi = result
            combined_barcode = f"{barcode1}{barcode2}{barcode3}"
            barcodes.append(combined_barcode)

            # Write the modified sequence to the output file
            y += 1
            outfile.write("@%s\n%s\n+\n%s\n" % (
                record.name + f"_{combined_barcode}_{umi}",
                seq,
                record.quality
            ))
        # Skip the read if result is None
        else:
            continue

outfile.close()

# Write out a list of whitelist barcodes
with open(args.whitelist, "w") as out_barcodes:
    for i in set(barcodes):
        out_barcodes.write("%s\n" % i)

# Log summary
log.write("The number of total reads: %s\n" % n)
log.write("The number of total reads that have a correct barcode and UMI: %s\n" % y)
log.write("The number of total recovered percent is: %s\n" % ((y / n) * 100))
log.close()
