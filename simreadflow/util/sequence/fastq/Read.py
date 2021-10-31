__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import os, sys
dis = os.path.abspath("../../../")
import gzip
sys.path.append(dis)
from Bio import SeqIO
from Path import to
import pyfastx


class read(object):

    def __init__(self):
        pass

    def fromgz(self, fastq_path, fastq_name, method='biopython'):
        names = []
        seqs = []
        placeholders = []
        qualities = []
        if method == 'biopython':
            with gzip.open(fastq_path + fastq_name + '.fastq.gz', 'rt') as handle:
                for record in SeqIO.parse(handle, 'fastq'):
                    # print()
                    seqs.append(''.join(record.seq))
                    names.append(''.join(record.name))
            return names, seqs, placeholders, qualities
        elif method == 'pyfastx':
            fq = pyfastx.Fastx(fastq_path + fastq_name + '.fastq.gz')
            for name, seq, qual, comment in fq:
                seqs.append(''.join(seq))
                names.append(''.join(name))
            return names, seqs, placeholders, qualities