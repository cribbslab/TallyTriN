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
L = logging.getLogger("complement_ployA.py")

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


tab = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


outfile = iotools.open_file(args.outname, "w")
log =  iotools.open_file(args.outname + ".log","w")
n = 0
y = 0
with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        y +=1
        if len(record.sequence) < 300:
            pass
        else:
            seq = record.sequence[50:300]
            m=regex.findall("(TTTTTTTTTTTTTTTTTTTT){e<=3}", str(seq))
            if m:
                n +=1
                sequence = reverse_complement_table(str(record.sequence))
                quality = str(record.quality)[::-1]
                outfile.write("@%s\n%s\n+\n%s\n" % (record.name, sequence, quality))
            else:
                seq = record.sequence[-30:-200:]
                m=regex.findall("(AAAAAAAAAAAAAAAAAAAA){e<=3}", str(seq))
                if m:
                    n +=1
                    outfile.write("@%s\n%s\n+\n%s\n" % (record.name, record.sequence, record.quality))
        

log.write("The number of total reads with polyA: %s\n" %(n))
log.write("The number of total reads is: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((n/y)*100))

log.close()
outfile.close()
