__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

from src.util.sequence.convert.ManyToSingle import manyToSingle as sm2s
from src.sequencing.reads.simulate.InitiatePool import initiatePool as simuip
from src.sequencing.reads.simulate.pcr.Amplify import amplify as pcr
from src.sequencing.reads.simulate.sequencing.Calling import calling as seq
from src.util.sequence.fastq.Write import write as wfastq
from Path import to


class starter(object):
    
    def __init__(self, *args, **kwargs):
        self.args = args[0]
        self.kwargs = kwargs
        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq

    def ondemand(self, ):
        # /*** block. Init a pool of sequences ***/
        init_seqs = simuip(
            cand_pool_fpn=self.args['init_seq_setting']['cand_pool_fpn'],
            cdna_fp=self.args['init_seq_setting']['cdna_fp'],
            cdna_num=self.args['init_seq_setting']['cdna_num'],
            umi_unit_pattern=self.args['init_seq_setting']['umi_unit_pattern'],
            umi_unit_len=self.args['init_seq_setting']['umi_unit_len'],
            is_seed=self.args['init_seq_setting']['is_seed'],
            prime_len=self.args['init_seq_setting']['prime_len'],
            umi_lib_fpn=self.args['init_seq_setting']['umi_lib_fpn'],
        ).generate()
        self.args['init_seqs'] = init_seqs
        print('Init pool of sequences has completed.')
        # print(init_seqs)
        # /*** block. PCR amplification ***/
        print('PCR amplification has started...')
        pcr = self.pcr(self.args).np()
        print('PCR amplification has completed.')
        # /*** block. Sequencing ***/
        print(pcr.keys())
        print('Sequencing has started...')
        pcr['seq_error'] = self.args['seq_error']
        seq = self.seq(pcr).np()
        print('Sequencing has completed.')
        print('Write seqs in fastq format')
        self.wfastq().togz(
            list_2d=seq['data'],
            sv_fp=self.args['write']['fastq_fp'],
            fn=self.args['write']['fastq_fn'],
        )
        return


if __name__ == "__main__":
    dis = '../../../'
    offset = '../' * 4 + dis
    DEFINE = {
        'fasta_fpn': 'data/omics/genomics/fasta/cdna/GRCh38/Homo_sapiens.GRCh38.cdna.all.fa',
        'cand_pool_fpn': offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt',
        'cdna_fp': offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
    }
    simu_params = {
        'init_seq_setting': {
            'cdna_num': 100,
            'umi_unit_pattern': 3,
            'umi_unit_len': 12,
            'is_seed': True,
            'prime_len': 20,
            'cand_pool_fpn': DEFINE['cand_pool_fpn'],
            'cdna_fp': DEFINE['cdna_fp'],
            'umi_lib_fpn': to('data/umi.txt'),
        },
        'ampl_rate': 0.85,
        'num_pcr': 10,
        'pcr_error': 1e-3,
        'seq_error': 1e-3,
        'use_seed': False,
        'seed': None,
        'write': {
            'fastq_fp': to('data/'),
            'fastq_fn': 'simu',
        }
    }
    p = starter(simu_params)
    print(p.ondemand())