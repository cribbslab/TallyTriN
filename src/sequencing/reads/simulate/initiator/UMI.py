__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import time
import os, sys
import numpy as np
dis = '../../../'
sys.path.append(os.path.abspath(dis))
from src.util.random.Sampling import sampling as ranspl
from src.util.random.Number import number as rannum
from src.util.file.read.Reader import reader as pfreader
from src.util.file.create.Folder import folder as crtfolder
from src.util.sequence.symbol.Single import single as dnasgl
from src.sequencing.reads.umi.Design import design as dumi
from src.sequencing.reads.similarity.distance.Hamming import hamming


class umi(object):

    def __init__(self, seq_num, is_seed=False, umi_unit_pattern=3, umi_unit_len=12, is_sv_umi_lib=True, umi_lib_fpn='./umi.txt', working_dir='./simu/', condis=['umi'], sim_thres=2, permutation=0):
        self.pfreader = pfreader()
        self.ranspl = ranspl()
        self.rannum = rannum()
        self.dnasgl = dnasgl()
        self.crtfolder = crtfolder()
        self.dumi = dumi
        self.is_seed = is_seed
        self.is_sv_umi_lib = is_sv_umi_lib
        self.umi_lib_fpn = umi_lib_fpn
        self.seq_num = seq_num
        self.umi_unit_pattern = umi_unit_pattern
        self.umi_unit_len = umi_unit_len
        self.condis = condis
        self.sim_thres = sim_thres
        self.permutation = permutation
        self.dna_map = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)
        self.crtfolder.osmkdir(working_dir)

    @property
    def umi_len(self, ):
        return self.umi_unit_pattern * self.umi_unit_len

    def pooling(self,):
        stime = time.time()
        seqs = []
        umi_pool = []
        umi_cnt = 0
        for id in np.arange(self.seq_num):
            read_struct_ref = {}
            if 'umi' in self.condis:
                umi_flag = False
                while not umi_flag:
                    umip = self.dumi(
                        dna_map=self.dna_map,
                        umi_unit_pattern=self.umi_unit_pattern,
                        pseudorandom_num=self.rannum.uniform(
                            low=0,
                            high=4,
                            num=self.umi_unit_len,
                            use_seed=self.is_seed,
                            seed=id + self.permutation * self.seq_num + umi_cnt,
                        ),
                    )
                    umi_i = umip.reoccur(is_sv=False)
                    edh = np.array([hamming().general(umi_i, j) for j in umi_pool])
                    # for j in umi_pool:
                    #     if hamming().general(umi_i, j) < self.sim_thres:
                    #         print(umi_i, j)
                    if len(edh[edh < self.sim_thres]) == 0:
                        # print(len(edh[edh < self.sim_thres]))
                        umi_pool.append(umi_i)
                        read_struct_ref['umi'] = umi_i
                        umi_flag = True
                        umip.write(res=umi_i, lib_fpn=self.umi_lib_fpn, is_sv=self.is_sv_umi_lib)
                    else:
                        # print(id)
                        umi_cnt += 1
            read_struct_pfd_order = {condi: read_struct_ref[condi] for condi in self.condis}
            seqs.append([self.paste([*read_struct_pfd_order.values()]), id, 'init'])
        print(umi_cnt)
        # print(umi_pool)

        etime = time.time()
        print("--->time for generating initial pool of sequences: {:.3f}s".format(etime-stime))
        return seqs

    def paste(self, read_struct=[]):
        return ''.join(read_struct)


if __name__ == "__main__":
    from Path import to
    DEFINE = {
        '': '',
    }
    # print(DEFINE['cand_pool_fpn'])
    p = umi(
        seq_num=1000,
        umi_unit_pattern=1,
        umi_unit_len=12,
        is_seed=True,

        is_sv_umi_lib=True,
        umi_lib_fpn=to('data/simu/umi/umi-plot.txt'),
        working_dir=to('data/simu/umi/'),

        condis=['umi'],
        sim_thres=3,
        permutation=1,
    )

    # print(p.umi_len)
    res = p.pooling()
    print(res)