__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import os, sys
import numpy as np
dis = '../../../'
sys.path.append(os.path.abspath(dis))
import time
from src.util.sequence.fasta.Read import read as sfasta
from src.util.random.Sampling import sampling as ranspl
from src.util.random.Number import number as rannum
from src.util.file.read.Reader import reader as pfreader
from src.util.file.create.Folder import folder as crtfolder
from src.util.sequence.symbol.Single import single as dnasgl
from src.sequencing.reads.umi.Design import design as umi
from src.sequencing.reads.seq.Design import design as seq


class general(object):

    def __init__(self, seq_num, is_seed=False, umi_unit_pattern=3, umi_unit_len=12, seq_len=20, is_sv_umi_lib=True, is_sv_seq_lib=True, umi_lib_fpn='./umi.txt', seq_lib_fpn='./seq.txt', working_dir='./simu/', permutation=0):
        self.pfreader = pfreader()
        self.ranspl = ranspl()
        self.rannum = rannum()
        self.dnasgl = dnasgl()
        self.sfasta = sfasta()
        self.crtfolder = crtfolder()
        self.umi = umi
        self.seq = seq
        self.is_seed = is_seed
        self.is_sv_umi_lib = is_sv_umi_lib
        self.is_sv_seq_lib = is_sv_seq_lib
        self.umi_lib_fpn = umi_lib_fpn
        self.seq_lib_fpn = seq_lib_fpn
        self.seq_num = seq_num
        self.umi_unit_pattern = umi_unit_pattern
        self.umi_unit_len = umi_unit_len
        self.seq_len = seq_len
        self.permutation = permutation
        self.dna_map = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)
        self.crtfolder.osmkdir(working_dir)

    @property
    def umi_len(self, ):
        return self.umi_unit_pattern * self.umi_unit_len

    def pooling(self, ):
        stime = time.time()
        seqs = [[self.paste(
            umi=self.umi(
                dna_map=self.dna_map,
                umi_unit_pattern=self.umi_unit_pattern,
                pseudorandom_num=self.rannum.uniform(
                    low=0,
                    high=4,
                    num=self.umi_unit_len,
                    use_seed=self.is_seed,
                    seed=id + self.permutation,
                ),
            ).general(lib_fpn=self.umi_lib_fpn, is_sv=self.is_sv_umi_lib),
            seq=self.seq(
                dna_map=self.dna_map,
                pseudorandom_num=self.rannum.uniform(
                    low=0,
                    high=4,
                    num=self.seq_len,
                    use_seed=self.is_seed,
                    seed=id + 50000 + self.permutation,
                ),
            ).general(lib_fpn=self.seq_lib_fpn, is_sv=self.is_sv_seq_lib),
        ), id, 'init'] for id in np.arange(self.seq_num)]
        etime = time.time()
        print("===>time for generating initial pool of sequences: {:.3f}s".format(etime-stime))
        return seqs

    def paste(self, umi, seq):
        return umi + seq


if __name__ == "__main__":
    DEFINE = {
        '': '',
    }
    # print(DEFINE['cand_pool_fpn'])
    p = general(
        seq_num=10,
        seq_len=20,
        umi_unit_pattern=1,
        umi_unit_len=12,
        is_seed=True,
    )

    # print(p.umi_len)
    res = p.pooling()
    print(res)