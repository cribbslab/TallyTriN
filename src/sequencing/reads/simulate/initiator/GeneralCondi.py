__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import os, sys
import numpy as np
dis = '../../../'
sys.path.append(os.path.abspath(dis))
from src.util.sequence.fasta.Read import read as sfasta
from src.util.random.Sampling import sampling as ranspl
from src.util.random.Number import number as rannum
from src.util.file.read.Reader import reader as pfreader
from src.util.file.create.Folder import folder as crtfolder
from src.util.sequence.symbol.Single import single as dnasgl
from src.sequencing.reads.umi.Design import design as umi
from src.sequencing.reads.seq.Design import design as seq
from src.sequencing.reads.spacer.Design import design as spacer
from src.sequencing.reads.primer.Design import design as primer


class generalCondi(object):

    def __init__(self, seq_num, is_seed=False, umi_unit_pattern=3, umi_unit_len=12, seq_len=20, spacer_len=4, primer_len=20, is_sv_umi_lib=True, is_sv_seq_lib=True, is_sv_spacer_lib=True, is_sv_primer_lib=True, umi_lib_fpn='./umi.txt', seq_lib_fpn='./seq.txt', spacer_lib_fpn='./spacer.txt', primer_lib_fpn='./primer.txt', working_dir='./simu/', condis=['umi', 'spacer', 'seq'], permutation=0):
        self.pfreader = pfreader()
        self.ranspl = ranspl()
        self.rannum = rannum()
        self.dnasgl = dnasgl()
        self.sfasta = sfasta()
        self.crtfolder = crtfolder()
        self.umi = umi
        self.seq = seq
        self.spacer = spacer
        self.primer = primer
        self.is_seed = is_seed
        self.is_sv_umi_lib = is_sv_umi_lib
        self.is_sv_seq_lib = is_sv_seq_lib
        self.is_sv_spacer_lib = is_sv_spacer_lib
        self.is_sv_primer_lib = is_sv_primer_lib
        self.umi_lib_fpn = umi_lib_fpn
        self.seq_lib_fpn = seq_lib_fpn
        self.spacer_lib_fpn = spacer_lib_fpn
        self.primer_lib_fpn = primer_lib_fpn
        self.seq_num = seq_num
        self.umi_unit_pattern = umi_unit_pattern
        self.umi_unit_len = umi_unit_len
        self.seq_len = seq_len
        self.spacer_len = spacer_len
        self.primer_len = primer_len
        self.condis = condis
        self.permutation = permutation
        self.dna_map = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)
        self.crtfolder.osmkdir(working_dir)

    @property
    def umi_len(self, ):
        return self.umi_unit_pattern * self.umi_unit_len

    def generate(self,):
        stime = time.time()
        seqs = []
        for id in np.arange(self.seq_num):
            read_struct_ref = {}
            if 'primer' in self.condis:
                primer_i = self.primer(
                    dna_map=self.dna_map,
                    pseudorandom_num=self.rannum.uniform(
                        low=0,
                        high=4,
                        num=self.primer_len,
                        use_seed=self.is_seed,
                        seed=id + 5000000*1 + self.permutation * self.seq_num,
                    ),
                ).general(lib_fpn=self.primer_lib_fpn, is_sv=self.is_sv_primer_lib)
                read_struct_ref['primer'] = primer_i
            if 'umi' in self.condis:
                umi_i = self.umi(
                    dna_map=self.dna_map,
                    umi_unit_pattern=self.umi_unit_pattern,
                    pseudorandom_num=self.rannum.uniform(
                        low=0,
                        high=4,
                        num=self.umi_unit_len,
                        use_seed=self.is_seed,
                        seed=id + self.permutation * self.seq_num,
                    ),
                ).general(lib_fpn=self.umi_lib_fpn, is_sv=self.is_sv_umi_lib)
                read_struct_ref['umi'] = umi_i
            if 'spacer' in self.condis:
                spacer_i = self.spacer(
                    dna_map=self.dna_map,
                    pseudorandom_num=self.rannum.uniform(
                        low=0,
                        high=4,
                        num=self.spacer_len,
                        use_seed=self.is_seed,
                        seed=id + 5000000*2 + self.permutation * self.seq_num,
                    ),
                ).general(lib_fpn=self.spacer_lib_fpn, is_sv=self.is_sv_spacer_lib)
                read_struct_ref['spacer'] = spacer_i
            if 'seq' in self.condis:
                seq_i = self.seq(
                    dna_map=self.dna_map,
                    pseudorandom_num=self.rannum.uniform(
                        low=0,
                        high=4,
                        num=self.seq_len,
                        use_seed=self.is_seed,
                        seed=id + 5000000*3 + self.permutation * self.seq_num,
                    ),
                ).general(lib_fpn=self.seq_lib_fpn, is_sv=self.is_sv_seq_lib)
                read_struct_ref['seq'] = seq_i
            read_struct_pfd_order = {condi: read_struct_ref[condi] for condi in self.condis}
            seqs.append([self.paste([*read_struct_pfd_order.values()]), id, 'init'])
        etime = time.time()
        print("--->time for generating initial pool of sequences: {:.3f}s".format(etime-stime))
        return seqs

    def paste(self, read_struct=[]):
        return ''.join(read_struct)


if __name__ == "__main__":
    DEFINE = {
        '': '',
    }
    # print(DEFINE['cand_pool_fpn'])
    p = generalCondi(
        seq_num=10,
        umi_unit_pattern=1,
        umi_unit_len=6,
        seq_len=20,
        spacer_len=4,
        primer_len=10,
        is_seed=True,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_spacer_lib=True,
        is_sv_primer_lib=True,
        umi_lib_fpn='./umi.txt',
        seq_lib_fpn='./seq.txt',
        spacer_lib_fpn='./spacer.txt',
        primer_lib_fpn='./primer.txt',
        working_dir='./simu/',

        condis=['umi', 'spacer', 'seq'],
        permutation=1,
    )

    # print(p.umi_len)
    res = p.generate()
    print(res)