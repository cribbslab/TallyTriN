import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import pandas as pd
import scipy.sparse as sparse
import scipy.io as io
import os


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("single_nucleotide_select.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='fastq file')
parser.add_argument("--outfile", default=None, type=str,
                    help='outfile fastq with trimers collapsed randomly')
args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #

outf =  iotools.open_file(args.outfile,"w")
log =  iotools.open_file(args.outfile + ".log","w")

with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:

        name = record.name.split("_")[0]
        barcode = record.name.split("_")[1]
        barcode = [barcode[i:i+3] for i in range(0, len(barcode), 3)]
        umi = record.name.split("_")[2]
        umi = [umi[i:i+3] for i in range(0, len(umi), 3)]

        barcode_collapse = []

        umi_collapse = []
        for trimer in barcode:
            barcode_collapse.append(list(trimer)[0])
            
        for trimer in umi:
            umi_collapse.append(list(trimer)[0])
        
        barcode_collapse = "".join(barcode_collapse)
        umi_collapse = "".join(umi_collapse)

        if len(barcode_collapse) == 10 and len(umi_collapse) == 6:

            outf.write("@%s\n%s\n+\n%s\n" % (name + "_" + barcode_collapse + "_" + umi_collapse, record.sequence, record.quality))
            

outf.close()

