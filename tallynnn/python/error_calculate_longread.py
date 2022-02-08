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
L = logging.getLogger("error_correct_longread.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--outfile", default=None, type=str,
                    help='output metrics file')

args = parser.parse_args()

L.info("args:")
print(args)

def most_common(lst):
    return max(set(lst), key=lst.count)


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
            #print(n)
            if n > 1:
                error_counter += 1
    return("".join(point_remove))


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
        umi = remove_point_mutations(umi)
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


# How many idels vs point mutations

umi_insertions = []
def allCharactersSame(s) :
    n = len(s)
    for i in range(1, n) :
        if s[i] != s[0] :
            return False
 
    return True

outfile = iotools.open_file(args.outfile, "w")

indel = 0
point = 0
total = 0
lines = 0
total_umis_counted = 0
with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        lines += 1
        seq_nano = record.sequence
        
        m=regex.finditer("(GTACTCTGCGTTGATACCACTGCTT){e<=0}", str(record.sequence))
        
        first = 0
        for i in m:
            first += 1
            if first == 1:
                total +=1
            umi = seq_nano[i.start()-30:i.start()]
            #print(umi)
            single_nuc = []

            trimers = [umi[i:i+3] for i in range(0, len(umi), 3)]

            n = 0
            error_counter = 0
            for trimer in trimers:

                if allCharactersSame(trimer):
                    pass
                if not allCharactersSame(trimer):
                    #print(trimer)
                    n += 1
                    #print(n)
                    if n >= 1:
                        error_counter += 1
            
            if error_counter ==0:
                total_umis_counted += 1
            if error_counter ==1:
                point += 1
                total_umis_counted += 1
                umi_insertions.append(umi)
            if error_counter ==2:
                total_umis_counted += 1
                point += 1
                umi_insertions.append(umi)
            if error_counter ==3:
                point += 1
                total_umis_counted += 1
                umi_insertions.append(umi)
        
            if error_counter > 3:
                indel += 1
                total_umis_counted += 1
                umi_insertions.append(umi)


outfile.write("Total UMIs that are free of errors: {}\n".format((lines-(point+indel))/lines *100))
outfile.write("Total UMIs that have either point or indels: {}\n".format((point+indel)/lines *100))
outfile.write("Total UMIs that have a point error: {}\n".format(point/lines *100))
outfile.write("Total UMIs that have either indels error: {}\n\n\n".format(indel/lines *100))
                        
        
corrected_umis = correct_umi(umi_insertions)

umis_not_correct = []

perfect = 0
total_umi = 0
for umi in corrected_umis:
    total_umi += 1
    if umi == 'GACTGCTACT':
        perfect += 1
    else:
        umis_not_correct.append(umi)

outfile.write("Of those that have errors, how many are then corrected as a percentage of total UMIs:\n\n")
outfile.write("Error UMIs that are free of errors: {}\n".format(perfect/total_umi *100))
outfile.write("Error UMIs that still have errors: {}\n".format((total_umi - perfect)/total_umi *100))
