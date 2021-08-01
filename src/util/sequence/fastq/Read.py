__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import os, sys
dis = os.path.abspath("../../../")
import gzip
print(dis)
sys.path.append(dis)
from Bio import SeqIO
from Path import to


class read(object):

    def __init__(self):
        pass

    def fromgz(self, fastq_path, fastq_name):
        seqs = []
        with gzip.open(fastq_path + fastq_name + '.fastq.gz', 'rt') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                # print()
                seqs.append(''.join(record.seq))
        # sequence = []
        # for seq in SeqIO.parse(fasta_path + fasta_name + '.fastq.gz', "fastq"):
        #     print(seq.seq)
        #     sequence.append(str(seq.seq))
        # sequence = ''.join(sequence)
        # if sequence == '':
        #     print('The sequence is empty.')
        return seqs
