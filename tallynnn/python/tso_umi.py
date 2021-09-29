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

outfile = iotools.open_file(args.outname, "w")
log =  iotools.open_file(args.outname + ".log","w")

n = 0
y = 0
with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        n += 1
        seq_nano = record.sequence
        
        m=regex.finditer("(AAGCAGTGGTATCAACGCAGAGTAAT){e<=3}", str(record.sequence))
        for i in m:
            after_smart = seq_nano[i.end()+1:]
            umi_tso = after_smart[:36]
            after_umi = after_smart[37:]
            record_new = record.name + "_" + str(umi_tso)
            quality_afterumismart = record.quality[i.end()+38:]
            seq_afterumismart = seq_nano[i.end()+38:]
            if len(umi_tso) == 36:
                y += 1
                # at the moment I am not including the umi for tso in the name
                outfile.write("@%s\n%s\n+\n%s\n" % (record_new, seq_afterumismart, quality_afterumismart))
            else:
                pass


log.write("The number of total reads: %s\n" %(n))
log.write("The number of total reads with a TSO UMI: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((y/n)*100))

log.close()
outfile.close()
