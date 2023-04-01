import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("correct_barcode_nano.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--outfile", default=None, type=str,
                    help='name for output fastq files')

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #


outf =  iotools.open_file(args.outfile,"w")
log =  iotools.open_file(args.outfile + ".log","w")


def most_common(lst):
    return max(set(lst), key=lst.count)


with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        
        seq = record.sequence
        
        name = record.name.split("_")[0]
        barcode = record.name.split("_")[1]
        barcode = [barcode[i:i+3] for i in range(0, len(barcode), 3)]
        umi = record.name.split("_")[2]
        umi = [umi[i:i+3] for i in range(0, len(umi), 3)]
            
        barcode_collapse = []
            
        umi_collapse = []
        for trimer in barcode:
            barcode_collapse.append(most_common(trimer))
            
        for trimer in umi:
            umi_collapse.append(most_common(trimer))
        
        barcode_collapse = "".join(barcode_collapse)
        barcode_collapse = barcode_collapse[:27]
        umi_collapse = "".join(umi_collapse)
        
        outf.write("@%s\n%s\n+\n%s\n" % (name + "_" + barcode_collapse + "_" + umi_collapse, record.sequence, record.quality))
            

outf.close()
