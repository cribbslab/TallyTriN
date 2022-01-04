import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import collections


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("polya_umi.py")

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

outfile = iotools.open_file(args.outname, "w")
log =  iotools.open_file(args.outname + ".log","w")

n = 0
y = 0
with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        n += 1
        seq_nano = record.sequence
        
        m=regex.finditer("(GTACTCTGCGTTGATACCACTGCTT){e<=0}", str(record.sequence))
        
        for i in m:
            after_polya = seq_nano[i.start()-18:]
            umi_polya = seq_nano[i.start()-18:i.start()]

            after_umi = seq_nano[:i.start()-18]

            record_new = record.name + str(umi_polya)
            
            quality_afterumipolya = record.quality[:i.start()-36]
            seq_afterumipolya = seq_nano[:i.start()-36]
            
            if len(umi_polya) == 18:
                y += 1
                outfile.write("@%s\n%s\n+\n%s\n" % (record_new, seq_afterumipolya, quality_afterumipolya))
            else:
                pass

log.write("The number of total reads: %s\n" %(n))
log.write("The number of total reads with a polyA UMI: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((y/n)*100))

log.close()
outfile.close()
