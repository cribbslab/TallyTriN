import sys
import logging
import glob

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("merge_counts.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #


infiles = glob.glob('*.counts.txt')
outfile = open('fileCounts.csv', 'w')

outfile.write("name,count\n")
for file in infiles:
    name = file.replace('.counts.txt','')
    file_open = open(file, 'r')
    for line in file_open:
        outfile.write("%s,%s" % (name, line))
outfile.close()
