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
parser.add_argument("--errors", default=None, type=str,
                    help='Number of errors to filter')
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


tab = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]

outfile = iotools.open_file("polya_revcomp.fastq", "w")
log =  iotools.open_file("revcomp" + ".log","w")
n = 0
y = 0
with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        y +=1
        if len(record.sequence) < 300:
            pass
        else:
            seq = record.sequence[30:300]
            m=regex.findall("(TTTTTTTTTTTTTTTTTTTT){e<=0}", str(seq))
            if m:
                n +=1
                sequence = reverse_complement_table(str(record.sequence))
                outfile.write("@%s\n%s\n+\n%s\n" % (record.name, sequence, record.quality))
            else:
                seq = record.sequence[-30:-200:]
                m=regex.findall("(AAAAAAAAAAAAAAAAAAAA){e<=0}", str(seq))
                if m:
                    n +=1
                    outfile.write("@%s\n%s\n+\n%s\n" % (record.name, record.sequence, record.quality))
                            
outfile.close()
log.close()

n = 0

collapsed_trimers = []
umis = []

with pysam.FastxFile("polya_revcomp.fastq") as fh:
    
    for record in fh:
        n += 1
        seq_nano = record.sequence
        
        m=regex.finditer("(GTACTCTGCGTTGATACCACTGCTT){e<=0}", str(record.sequence))
        if n < 100000:
            for i in m:
                umi = seq_nano[i.start()-30:i.start()]
                #print(umi)
                single_nuc = []
                trimers = [umi[i:i+3] for i in range(0, len(umi), 3)]

                for trimer in trimers:
                    single_nuc.append(trimer)

                collapse = "".join(single_nuc)
                collapsed_trimers.append(collapse)
                umis.append(umi)
        else:
            break

outfile = iotools.open_file(args.outfile, "w")
    
# count the number of perfect 
from collections import Counter
counts = []
n = 0
for sequence in collapsed_trimers:
    n += 1
    count = sum(1 for a, b in zip("GGGAAACCCTTTGGGCCCTTTAAACCCTTT", sequence) if a != b)
    #print(count)
    counts.append(count)

outfile.write("The percent of perfect UMIs with trimers: %s\n"%(Counter(counts)[0]/n*100))
outfile.write("The percent of one missmatch UMI with trimers: %s\n"%(Counter(counts)[1]/n*100))
outfile.write("The percent of two missmatch UMI with trimers: %s\n"%(Counter(counts)[2]/n*100))
outfile.write("The percent of three missmatch UMI with trimers: %s\n"%(Counter(counts)[3]/n*100))
outfile.write("The percent of 4 missmatch UMI with trimers: %s\n"%(Counter(counts)[4]/n*100))
outfile.write("The percent of 5 missmatch UMI with trimers: %s\n"%(Counter(counts)[5]/n*100))
outfile.write("The percent of 6 missmatch UMI with trimers: %s\n"%(Counter(counts)[6]/n*100))
outfile.write("The percent of 7 missmatch UMI with trimers: %s\n"%(Counter(counts)[7]/n*100))
outfile.write("The percent of 8 missmatch UMI with trimers: %s\n"%(Counter(counts)[8]/n*100))
outfile.write("The percent of 9 missmatch UMI with trimers: %s\n"%(Counter(counts)[9]/n*100))
outfile.write("The percent of 10 missmatch UMI with trimers: %s\n\n"%(Counter(counts)[10]/n*100))


# Calculate the recovery of each umi that isnt perfect

umi_insertions = []
for collapse, umi in zip(collapsed_trimers, umis):
    umi_insertions.append(umi)


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
        num_errors = 0
        new_umi = []
        umi, errors = remove_point_mutations(umi)
        if errors > int(args.errors):
            num_errors =+ 1
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
    return(corrected_umis, num_errors)


corrected_umis, num_errors = correct_umi(umi_insertions)



umis_not_correct = []

perfect = 0
total = 0
for umi in corrected_umis:
    total += 1
    if umi == 'GACTGCTACT':
        perfect += 1
    else:
        umis_not_correct.append(umi)


outfile.write("Total UMIs: {}\n".format(total))
outfile.write("Total perfect UMIs: {}\n".format(perfect))
outfile.write("Number of UMIs that have errors greater than {}: {}\n".format(args.errors, num_errors))
outfile.write("Total UMIs that are free of errors: {}\n".format(perfect/total *100))
outfile.write("Total UMIs that have errors: {}\n".format((total - perfect)/total *100))
