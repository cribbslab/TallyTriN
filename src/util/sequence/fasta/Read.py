__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import os, sys
dis = os.path.abspath("../../../../")
print(dis)
sys.path.append(dis)
from Bio import SeqIO
from src.sequencing.reads.umi.Library import library as liblogginger


class read(object):

    def __init__(self):
        pass

    @liblogginger()
    def seqIO(self, fasta_path, fasta_name, lib_fpn='./seq.txt', is_sv=True):
        sequence = []
        for seq in SeqIO.parse(fasta_path + fasta_name + '.fasta', "fasta"):
            # print(seq.seq)
            sequence.append(str(seq.seq))
        sequence = ''.join(sequence)
        if sequence == '':
            print('The sequence is empty.')
        return sequence

    def save_(self, list_2d, sv_fp):
        for i, e in enumerate(list_2d):
            prot_name = str(e[0])
            seq = str(e[1])
            print('No.{} saving {} in FASTA format.'.format(i+1, prot_name))
            f = open(sv_fp + prot_name + '.fasta', 'w')
            f.write('>' + prot_name + '\n')
            f.write(seq + '\n')
            f.close()
        return 0


if __name__ == "__main__":
    p = read()
    from Path import to
    offset = '../' * 7
    print(offset)
    print(p.get(
        fasta_path=offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
        fasta_name='ENST00000379435.3',
    ))