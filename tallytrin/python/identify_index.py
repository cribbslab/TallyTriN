import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
from collections import Counter
import collections
import os.path
import gzip

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("identify_index.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--name", default=None, type=str,
                    help='name of sample')

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


name = args.name


dict = {'AAATTTGGGCCC': '1',
        'TTTCCCAAAGGG': '2',
        'GGGAAACCCTTT': '3',
        'CCCGGGTTTAAA': '4',
        'AAACCCGGGAAA': '5',
        'TTTGGGAAATTT': '6',
        'GGGTTTCCCGGG': '7',
        'CCCAAATTTCCC': '8',
        'AAAGGGAAAGGG': '9',
        'TTTAAATTTAAA': '10',
        'GGGCCCGGGCCC': '11',
        'CCCTTTCCCTTT': '12'}




with pysam.FastxFile(args.infile) as fh:
    
    for record in fh:
        
        seq_nano = record.sequence

        m=regex.finditer("(AAGCAGTGGT){e<=1}", str(seq_nano))
        
        
        for i in m:
            barcode = seq_nano[int(i.start())-12:int(i.start())]
            barcode = remove_point_mutations(barcode)[0]
            barcode = barcode[::3]
            if barcode == 'ATGC':
                fname = name + '' + "Sample1.fastq"
            elif barcode == 'TCAG':
                fname = name + '' + "Sample2.fastq"
            elif barcode == 'GACT':
                fname = name + '' + "Sample3.fastq"
            elif barcode == 'CGTA':
                fname = name + '' + "Sample4.fastq"
            elif barcode == 'ACGA':
                fname = name + '' + "Sample5.fastq"
            elif barcode == 'TGAT':
                fname = name + '' + "Sample6.fastq"
            elif barcode == 'GTCG':
                fname = name + '' + "Sample7.fastq"
            elif barcode == 'CATC':
                fname = name + '' + "Sample8.fastq"
            elif barcode == 'AGAG':
                fname = name + '' + "Sample9.fastq"
            elif barcode == 'TATA':
                fname = name + '' + "Sample10.fastq"
            elif barcode == 'GCGC':
                fname = name + '' + "Sample11.fastq"
            elif barcode == 'CTCT':
                fname = name + '' + "Sample12.fastq"
            else:
                fname = name + '' + "Unidentified.fastq"
            
            if fname == name + '' + 'Unidentified.fastq':
                m=regex.finditer("(ACCACTGCTT){e<=1}", str(seq_nano))
                for i in m:
                    barcode = seq_nano[int(i.end()):int(i.end()+12)]
                    barcode = remove_point_mutations(barcode)[0]
                    barcode = barcode[::3]
                    if barcode == 'GCAT':
                        fname = name + '' + "Sample1.fastq"
                    elif barcode == 'CTGA':
                        fname = name + '' + "Sample2.fastq"
                    elif barcode == 'AGTC':
                        fname = name + '' + "Sample3.fastq"
                    elif barcode == 'TACG':
                        fname = name + '' + "Sample4.fastq"
                    elif barcode == 'TCGT':
                        fname = name + '' + "Sample5.fastq"
                    elif barcode == 'ATCA':
                        fname = name + '' + "Sample6.fastq"
                    elif barcode == 'CGAC':
                        fname = name + '' + "Sample7.fastq"
                    elif barcode == 'GATG':
                        fname = name + '' + "Sample8.fastq"
                    elif barcode == 'CTCT':
                        fname = name + '' + "Sample9.fastq"
                    elif barcode == 'TATA':
                        fname = name + '' + "Sample10.fastq"
                    elif barcode == 'GCGC':
                        fname = name + '' + "Sample11.fastq"
                    elif barcode == 'CTCT':
                        fname = name + '' + "Sample12.fastq"
                    else:
                        fname = name + '' + "Unidentified.fastq"
                if os.path.exists('seperate_samples.dir/' + fname):

                    with open('seperate_samples.dir/'+ fname, "a") as myfile:
                        myfile.write("@%s\n%s\n+\n%s\n" % (record.name, record.sequence, record.quality))
                else:
                    with iotools.open_file('seperate_samples.dir/' + fname, "w") as myfile:
                        myfile.write("@%s\n%s\n+\n%s\n" % (record.name, record.sequence, record.quality))

            else:
                break


