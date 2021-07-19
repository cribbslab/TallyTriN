__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import os, sys
dis = '../../'
sys.path.append(os.path.abspath(dis))
import time
from numba import jit
from src.util.sequence.Fasta import fasta as sfasta
from src.util.sampling.Group import group as spgroup
from src.util.sampling.Single import single as spsingle
from src.util.file.read.Reader import reader as pfreader
from src.util.sequence.symbol.Single import single as dnasgl
from Path import to


class native(object):

    def __init__(self, cand_pool_fpn, cdna_fp, cdna_num, is_seed_cdna=False, umi_pattern=3, umi_unit_len=12, prime_len=20):
        self.pfreader = pfreader()
        self.spgroup = spgroup()
        self.spsingle = spsingle()
        self.dnasgl = dnasgl()
        self.sfasta = sfasta()
        self.is_seed_cdna = is_seed_cdna
        self.cdna_fp = cdna_fp
        self.cand_pool_fpn = cand_pool_fpn
        self.df_cand_pool = self.pfreader.generic(df_fpn=self.cand_pool_fpn)
        self.cdna_num = cdna_num
        self.umi_pattern = umi_pattern
        self.umi_unit_len = umi_unit_len
        self.prime_len = prime_len
        self.umi_len = self.umi_pattern * self.umi_unit_len
        self.dna_map = self.dnasgl.todict(bases=self.dnasgl.get(universal=True), reverse=True)

    @jit()
    def generate(self, ):
        stime = time.time()
        df_cand_sel = self.spgroup.uniform(
            df=self.df_cand_pool,
            num=self.cdna_num,
            use_seed=self.is_seed_cdna,
            seed=1,
            replace=True,
        )
        seqs = [self.paste(
            seq=self.sfasta.get(
                fasta_path=self.cdna_fp,
                fasta_name=cand
            ),
            umi=self.umiSimu(),
            prime=self.primerSimu(),
        ) for cand in df_cand_sel.values.squeeze()]
        etime = time.time()
        print("time: {}s".format(etime-stime))
        return seqs

    def umiSimu(self, ):
        atcg_simu = ''.join([
            self.dna_map[i]*self.umi_pattern for i in self.spsingle.uniform(low=0, high=4, num=self.umi_unit_len)
        ])
        return atcg_simu

    def primerSimu(self, ):
        atcg_simu = ''.join([
            self.dna_map[i] for i in self.spsingle.uniform(low=0, high=4, num=self.prime_len)
        ])
        return atcg_simu

    def paste(self, seq, umi, prime):
        return prime + umi + seq + umi + ''.join(list(reversed(prime)))


if __name__ == "__main__":
    offset = '../' * 3 + dis
    DEFINE = {
        'cand_pool_fpn': offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt',
        'cdna_fp': offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
    }
    p = native(
        cand_pool_fpn=DEFINE['cand_pool_fpn'],
        cdna_fp=DEFINE['cdna_fp'],
        cdna_num=10
    )
    p.generate()
