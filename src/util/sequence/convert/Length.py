__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__author__ = "Adam Cribbs lab"

import re
import sys

import pandas as pd

sys.path.append('../../../../')
from Bio import SeqIO
from src.util.file.write.Writer import writer as pfwriter
from src.util.file.read.Reader import reader as pfreader
from src.util.sequence.fasta.Read import read as sfasta
from Path import to


class length:

    def __init__(self, list_fpn=None, ):
        self.pfreader = pfreader()
        self.pfwriter = pfwriter()
        self.prot_df = self.pfreader.generic(list_fpn)
        print(self.prot_df)

    def sequence(self, fasta_path, thres=150, sv_fpn=None, is_sv=False):
        t = []
        for i in self.prot_df.index:
            prot_name = self.prot_df.iloc[i, 0]
            if i % 100 == 0:
                print('No.{} protein {}'.format(i, prot_name))
            sequence = sfasta().seqIO(
                fasta_path=fasta_path,
                fasta_name=prot_name,
                is_sv=False,
                lib_fpn=None
            )
            # print(sequence)
            if len(sequence) < thres:
                t.append(prot_name)
        print(len(t))
        if is_sv:
            self.pfwriter.generic(pd.DataFrame(t), sv_fpn=sv_fpn)
        return self.prot_df


if __name__ == "__main__":

    offset = '../' * 7
    DEFINE = {
        'cand_pool_fpn': offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt',
        'cdna_fp': offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
    }
    p = length(list_fpn=DEFINE['cand_pool_fpn'])
    print(p.sequence(
        fasta_path=DEFINE['cdna_fp'],
        thres=150,
        is_sv=True,
        sv_fpn=offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n_l150.txt',
    ))