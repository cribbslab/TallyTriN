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
        # eded = []
        # for i, u1 in enumerate(umi_pool[:50]):
        #     for j, u2 in enumerate(umi_pool[:50]):
        #         eded.append([i, j, hamming().general(u1, u2)])
        # import pandas as pd
        # import seaborn as sns
        # import matplotlib.pyplot as plt
        # tr = pd.DataFrame(eded, columns=['UMI1', 'UMI2', 'Edit distance'])
        #
        # ax = sns.kdeplot(data=tr.loc[tr['Edit distance'] > 0], x='Edit distance', linewidth=4, fill=True, common_norm=False, palette="crest", alpha=0.4,)
        # from matplotlib import rcParams
        # rcParams['font.family'] = 'sans-serif'
        # rcParams['font.sans-serif'] = ['Tahoma']
        # ax.set_xlabel('Edit distance', fontsize=13)
        # ax.set_ylabel('Density', fontsize=13)
        # ax.spines['right'].set_color('none')
        # ax.spines['top'].set_color('none')
        # plt.show()
        #
        # ty = tr['Edit distance'].max()
        # tr['Corr'] = tr['Edit distance'].apply(lambda x: 1-x/ty)
        # print(tr)
        #
        # # print(tr.loc[:100, :])
        # g = sns.relplot(
        #     data=tr,
        #     x='UMI1', y='UMI2', hue='Corr', size='Corr',
        #     palette="vlag", hue_norm=(-1, 1), edgecolor=".7",
        #     height=10, sizes=(50/1.5, 250/1.5), size_norm=(-.2, .8),
        # )
        #
        # # Tweak the figure to finalize
        # g.set(xlabel="", ylabel="", aspect="equal")
        # g.despine(left=True, bottom=True)
        # g.ax.margins(.02)
        # for label in g.ax.get_xticklabels():
        #     label.set_rotation(90)
        # for artist in g.legend.legendHandles:
        #     artist.set_edgecolor(".7")
        # plt.show()


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