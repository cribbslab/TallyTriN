import sys
import regex
import cgatcore.iotools as iotools
import pysam
from collections import Counter
import pandas as pd
import logging
import argparse
import re



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("trans_count.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='input bam  file')
parser.add_argument("--outfile", default=None, type=str,
                    help='name for output counts file')

args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #


bamfile = pysam.AlignmentFile(args.infile)

trans_list = []


for line in bamfile:
    transcript = line.reference_name
    if transcript is not None:
        trans_list.append(transcript)
    
counter_trans = Counter(trans_list)    

df = pd.DataFrame.from_dict(counter_trans, orient='index').reset_index()

df.set_index("index", inplace=True)

df.to_csv(args.outfile, sep="\t")

bamfile.close() 
