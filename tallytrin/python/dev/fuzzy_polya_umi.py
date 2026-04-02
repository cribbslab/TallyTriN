import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import collections
import parasail


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
parser.add_argument("--errors", default=None, type=str,
                    help='Number of errors to remove reads')
args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #



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
            n +=1
        if not allCharactersSame(trimer):
            base = collections.Counter(trimer).most_common(1)[0][0]
            base = base+base+base
            point_remove.append(base)
            error_counter += 1
    return("".join(point_remove), error_counter)


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


outfile = iotools.open_file(args.outname, "w")
log =  iotools.open_file(args.outname + ".log","w")

n = 0
y = 0

# ------------------ Adapter Pattern ------------------
adapter_pattern = "GTACTCTGCGTTGATACCACTGCTT"

# ------------------ Substitution Matrix ------------------
# ACGT substitution matrix: 2 for match, -1 for mismatch, 1 for gap penalty
matrix = parasail.matrix_create("ACGT", 2, -1)

# ------------------ Score Threshold ------------------
# Formula: (1 - error) * primer_length * match_score
error = 0.35
primer_length = len(adapter_pattern)  # 25
match_score = 2
min_score_threshold = (1 - error) * primer_length * match_score  # 32.5


with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        n += 1
        seq_nano = record.sequence[-200:]
        result = parasail.sw_trace_striped_16(adapter_pattern, seq_nano, 1, 1, matrix)  # Gap open/extend penalties
        #m=regex.finditer("(GTACTCTGCGTTGATACCACTGCTT){e<=0}", seq_nano)
        if result is None or result.score is None or result.traceback is None:
                    continue  # Skip this pattern if alignment failed
                  
        if result.score >= min_score_threshold:
            aligned_length = len(result.traceback.ref.replace("-", ""))
            match_pos = result.end_ref - aligned_length + 1  # Start position of the match
            m = [match_pos]

        for i in m:

            umi_polya = seq_nano[i-30:i]
            umi_polya = seq_nano[i-36:i]
            umi_polya, errors = remove_point_mutations(umi_polya)
            
            if errors > int(args.errors):
                pass
            else:


                record_new = record.name + "_" + str(umi_polya)
            
                #if len(umi_polya) == 30:
                if len(umi_polya) == 36:
                    y += 1
                    outfile.write("@%s\n%s\n+\n%s\n" % (record_new, record.sequence, record.quality))
                else:
                    pass

log.write("The number of total reads: %s\n" %(n))
log.write("The number of total reads with a polyA UMI: %s\n" %(y))
log.write("The number of total recovered percent is: %s\n" %((y/n)*100))

log.close()
outfile.close()
