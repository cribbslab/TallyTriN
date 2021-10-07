__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import os, sys
dis = '../../../'
sys.path.append(os.path.abspath(dis))
import time
# from numba import jit
from src.util.sequence.fasta.Read import read as sfasta
from src.util.random.Sampling import sampling as ranspl
from src.util.random.Number import number as rannum
from src.util.file.read.Reader import reader as pfreader
from src.util.file.create.Folder import folder as crtfolder
from src.util.sequence.symbol.Single import single as dnasgl
from src.sequencing.reads.umi.Design import design as umi
from src.sequencing.reads.primer.Design import design as primer


class cdna(object):

    def __init__(self, cand_pool_fpn, cdna_fp, cdna_num, is_seed=False,is_sv_umi_lib=True, is_sv_primer_lib=True, is_sv_seq_lib=True, umi_unit_pattern=3, umi_unit_len=12, primer_len=20, seq_lib_fpn='./seq.txt', umi_lib_fpn='./umi.txt', primer_lib_fpn='./primer.txt', working_dir='./simu/'):
        self.pfreader = pfreader()
        self.ranspl = ranspl()
        self.rannum = rannum()
        self.dnasgl = dnasgl()
        self.sfasta = sfasta()
        self.crtfolder = crtfolder()
        self.umi = umi
        self.primer = primer
        self.is_sv_umi_lib = is_sv_umi_lib
        self.is_sv_primer_lib = is_sv_primer_lib
        self.is_sv_seq_lib = is_sv_seq_lib
        self.is_seed = is_seed
        self.cdna_fp = cdna_fp
        self.cand_pool_fpn = cand_pool_fpn
        self.seq_lib_fpn = seq_lib_fpn
        self.umi_lib_fpn = umi_lib_fpn
        self.primer_lib_fpn = primer_lib_fpn
        self.df_cand_pool = self.pfreader.generic(df_fpn=self.cand_pool_fpn)
        self.cdna_num = cdna_num
        self.umi_unit_pattern = umi_unit_pattern
        self.umi_unit_len = umi_unit_len
        self.primer_len = primer_len
        self.dna_map = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)
        self.crtfolder.osmkdir(working_dir)

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
            seq=self.sfasta.seqIO(
                fasta_path=self.cdna_fp,
                fasta_name=cand,
                is_sv=self.is_sv_seq_lib,
                lib_fpn=self.seq_lib_fpn
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
            ).reoccur(lib_fpn=self.umi_lib_fpn, is_sv=self.is_sv_umi_lib),
            primer=self.primer(
                dna_map=self.dna_map,
                pseudorandom_num=self.rannum.uniform(
                    low=0,
                    high=4,
                    num=self.primer_len,
                    use_seed=self.is_seed,
                    seed=id,
                ),
            # ).general(primer_lib_fpn=self.primer_lib_fpn, is_sv=self.is_sv_primer_lib),
            ).tsoatdbio(lib_fpn=self.primer_lib_fpn, is_sv=self.is_sv_primer_lib),
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
    p = cdna(
        cand_pool_fpn=DEFINE['cand_pool_fpn'],
        cdna_fp=DEFINE['cdna_fp'],
        cdna_num=10,
        umi_unit_pattern=3,
        umi_unit_len=12,
        is_seed=True,
        primer_len=20,
    )

    # print(p.umi_len)
    res = p.generate()
    print(res)