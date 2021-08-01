__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import os, sys
dis = '../../../'
sys.path.append(os.path.abspath(dis))
import time
# from numba import jit
from src.util.sequence.Fasta import fasta as sfasta
from src.util.random.Sampling import sampling as ranspl
from src.util.random.Number import number as rannum
from src.util.file.read.Reader import reader as pfreader
from src.util.sequence.symbol.Single import single as dnasgl
from src.sequencing.reads.umi.Design import design as umi
from src.sequencing.reads.primer.Design import design as primer


class initiatePool(object):

    def __init__(self, cand_pool_fpn, cdna_fp, cdna_num, is_seed=False, umi_unit_pattern=3, umi_unit_len=12, prime_len=20, umi_lib_fpn='./umi.txt', primer_lib_fpn='./primer.txt'):
        self.pfreader = pfreader()
        self.ranspl = ranspl()
        self.rannum = rannum()
        self.dnasgl = dnasgl()
        self.sfasta = sfasta()
        self.umi = umi
        self.primer = primer
        self.is_seed = is_seed
        self.cdna_fp = cdna_fp
        self.cand_pool_fpn = cand_pool_fpn
        self.umi_lib_fpn = umi_lib_fpn
        self.primer_lib_fpn = primer_lib_fpn
        self.df_cand_pool = self.pfreader.generic(df_fpn=self.cand_pool_fpn)
        self.cdna_num = cdna_num
        self.umi_unit_pattern = umi_unit_pattern
        self.umi_unit_len = umi_unit_len
        self.prime_len = prime_len
        self.dna_map = self.dnasgl.todict(bases=self.dnasgl.get(universal=True), reverse=True)

    @property
    def umi_len(self, ):
        return self.umi_unit_pattern * self.umi_unit_len

    def generate(self, ):
        stime = time.time()
        df_cand_sel = self.ranspl.uniform(
            data=self.df_cand_pool,
            num=self.cdna_num,
            use_seed=self.is_seed,
            seed=1,
            replace=True,
        )
        # print(df_cand_sel)
        seqs = [[self.paste(
            seq=self.sfasta.get(
                fasta_path=self.cdna_fp,
                fasta_name=cand,
            ),
            umi=self.umi(
                dna_map=self.dna_map,
                umi_unit_pattern=self.umi_unit_pattern,
                pseudorandom_num=self.rannum.uniform(
                    low=0,
                    high=4,
                    num=self.umi_unit_len,
                    use_seed=self.is_seed,
                    seed=id,
                ),
            ).general(umi_lib_fpn=self.umi_lib_fpn),
            primer=self.primer(
                dna_map=self.dna_map,
                pseudorandom_num=self.rannum.uniform(
                    low=0,
                    high=4,
                    num=self.prime_len,
                    use_seed=self.is_seed,
                    seed=id,
                ),
            ).general(primer_lib_fpn=self.primer_lib_fpn),
        ), id, 'init'] for id, cand in enumerate(df_cand_sel.values.squeeze())]
        etime = time.time()
        print("time: {}s".format(etime-stime))
        return seqs

    def paste(self, seq, umi, primer):
        return primer + umi + seq + umi + ''.join(list(reversed(primer)))


if __name__ == "__main__":
    offset = '../' * 5 + dis
    DEFINE = {
        'cand_pool_fpn': offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt',
        'cdna_fp': offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
    }
    # print(DEFINE['cand_pool_fpn'])
    p = initiatePool(
        cand_pool_fpn=DEFINE['cand_pool_fpn'],
        cdna_fp=DEFINE['cdna_fp'],
        cdna_num=10,
        umi_unit_pattern=3,
        umi_unit_len=12,
        is_seed=True,
        prime_len=20,
    )

    # print(p.umi_len)
    res = p.generate()
    print(res)