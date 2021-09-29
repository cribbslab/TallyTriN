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
L = logging.getLogger("tso_umi.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')

args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #

log =  iotools.open_file(args.outname + ".log","w")

bamfile = pysam.AlignmentFile(args.infile, "rb")
split = pysam.AlignmentFile(args.outname, "wb", template=bamfile)

for line in bamfile:
    if line.has_tag("SA"):
        split.write(line)
    else:
        pass

bamfile.close()
split.close()


log.close()
