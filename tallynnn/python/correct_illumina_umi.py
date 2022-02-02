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
L = logging.getLogger("correct_illumina_umi.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--read1", default=None, type=str,
                    help='infile fastq.1.gz  file')
parser.add_argument("--read2", default=None, type=str,
                    help='infile fastq.2.gz  file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')

args = parser.parse_args()

L.info("args:")
print(args)

log =  iotools.open_file(args.outname + ".log","w")

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #

import collections

def allCharactersSame(s) :
    n = len(s)
    for i in range(1, n) :
        if s[i] != s[0] :
            return False
 
    return True

def remove_point_mutations(umi):
    
    point_remove = []
    
    trimers = [umi[i:i+3] for i in range(0, len(umi), 3)]

    n = 0
    error_counter = 0
    for trimer in trimers:

        if allCharactersSame(trimer):
            point_remove.append(trimer)
        if not allCharactersSame(trimer):
            base = collections.Counter(trimer).most_common(1)[0][0]
            base = base+base+base
            point_remove.append(base)
            n += 1
    return("".join(point_remove), n)


def remove_indels(x, umi, first):
    
    substr = umi[x:x+3]
    
    if ("CCC" in substr or "GGG" in substr or "AAA" in substr or "TTT" in substr) and first:
        pass
    
#    else:
#        # correct for point mutations at beginning
#        base = collections.Counter(substr).most_common(1)[0][0]
#        substr = base+base+base
#        first = False
        
    
    if not ("CCC" in substr or "GGG" in substr or "AAA" in substr or "TTT" in substr):
        # if there is an error then check to see if the next base is a perfect match
        substr = umi[x+1:x+4]
        
        if ("CCC" in substr or "GGG" in substr or "AAA" in substr or "TTT" in substr):
            # then check to see if perfect trimer
            pass

        elif not ("CCC" in substr or "GGG" in substr or "AAA" in substr or "TTT" in substr):
            # if not then check to see if there is a 
            substr = umi[x+2:x+5]

            
            
    return(substr)




def correct_umi(umis):
    n = 0
    corrected_umis = []
    for umi in umis:
        n +=1

        new_umi = []
        umi, errors = remove_point_mutations(umi)
        if errors > 2:
            pass
        else:
            for x in range(0, len(umi)):
                if x % 3 == 0:

                    if x == 0:
                        #print(x, umi)
                        sub_umi = remove_indels(x, umi, first=True)
                    else:
                        sub_umi = remove_indels(x, umi, first=False)

                    new_umi.append(sub_umi)

            new_umi = "".join(new_umi)
            final_umi = new_umi[:30]
            corrected_umis.append(final_umi[::3])
    return(corrected_umis)


outfile = iotools.open_file(args.outname, "w")
n = 0
y = 0
with pysam.FastxFile(args.read1) as fh, pysam.FastxFile(args.read2) as fh2:
    
    for record, record2 in zip(fh, fh2):
        n +=1
        umi = record.sequence[0:30]
        #print(umi)
        if "N" in umi:
            pass
        else:
            y += 1
            single_nuc = []
            trimers = [umi[i:i+3] for i in range(0, len(umi), 3)]
            collapsed_trimer = "".join(correct_umi(trimers))
            
            record_new = record2.name + "_" + str(collapsed_trimer)
            
            outfile.write("@%s\n%s\n+\n%s\n" % (record_new, record2.sequence, record2.quality))

outfile.close()


log.write("The number of total reads: %s\n" %(n))
log.write("The number of total reads with a corrected UMI: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((y/n)*100))

log.close()
outfile.close()
