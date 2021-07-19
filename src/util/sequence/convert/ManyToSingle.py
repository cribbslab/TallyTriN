__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__author__ = "Adam Cribbs lab"

import re
import sys
sys.path.append('../../../../')
from Bio import SeqIO
from src.util.file.write.Writer import writer as pfwriter
from Path import to


class manyToSingle(object):

    def __init__(self, fasta_fpn):
        self.pfwriter = pfwriter()
        self.fasta_fpn = fasta_fpn
        self.fas_id = []
        self.fas_seq = []
        self.fas_name = []
        self.fas_dpt = []
        for fas in SeqIO.parse(self.fasta_fpn, "fasta"):
            self.fas_id.append(fas.id)
            self.fas_seq.append(fas.seq)
            self.fas_name.append(fas.name)
            self.fas_dpt.append(fas.description)
        # print(self.fas_seq)

    def filter1(self, fasta_name):
        fasta_name_ = re.split('\|', fasta_name)[1]
        return fasta_name_

    def save(self, sv_fp):
        target_ids = self.getTargetId()
        print(target_ids)
        target_len = len(target_ids)
        for i in range(target_len):
            f = open(sv_fp + target_ids[i] + '.fasta', 'w')
            f.write('>' + target_ids[i] + '\n')
            f.write(str(self.fas_seq[i]))
            f.close()
        return

    def svid(self, sv_fpn):
        return self.pfwriter.generic(self.fas_id, sv_fpn=sv_fpn)

    def getTargetId(self, ):
        target_ids = []
        for _, id in enumerate(self.fas_id):
            target_ids.append(re.sub(r'^.*\|', "", str(id)))
            # target_ids.append(self.filter1(id))
        return target_ids


if __name__ == "__main__":
    DEFINE = {
        'normal':{
            'fasta_fpn': to('data/omics/genomics/fasta/cdna/GRCh38/Homo_sapiens.GRCh38.cdna.all.fa'),
            'sv_fpn_ids': to('data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt'),
            'sv_fasta_fp': to('data/omics/genomics/fasta/cdna/GRCh38/'),
        },
    }

    p = manyToSingle(
        fasta_fpn=DEFINE['normal']['fasta_fpn']
    )

    # print(p.fas_name)
    # print(len(p.fas_name))

    # print(p.svid(sv_fpn=DEFINE['normal']['sv_fpn_ids']))

    # print(p.save(sv_fp=DEFINE['normal']['sv_fasta_fp']))