__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from src.sequencing.reads.simulate.initiator.CDNA import cdna as simuip
from src.sequencing.reads.simulate.pcr.Amplify import amplify as pcr
from src.sequencing.reads.simulate.sequencing.Calling import calling as seq
from src.util.sequence.fastq.Write import write as wfastq
from Path import to


class cdna(object):
    
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
            primer_len=self.args['init_seq_setting']['primer_len'],
            is_sv_seq_lib=self.args['init_seq_setting']['is_sv_seq_lib'],
            is_sv_umi_lib=self.args['init_seq_setting']['is_sv_umi_lib'],
            is_sv_primer_lib=self.args['init_seq_setting']['is_sv_primer_lib'],
            seq_lib_fpn=self.args['init_seq_setting']['seq_lib_fpn'],
            umi_lib_fpn=self.args['init_seq_setting']['umi_lib_fpn'],
            primer_lib_fpn=self.args['init_seq_setting']['primer_lib_fpn'],
            working_dir=self.args['init_seq_setting']['working_dir'],
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
    offset = '../' * 5 + dis
    DEFINE = {
        'fasta_fpn': 'data/omics/genomics/fasta/cdna/GRCh38/Homo_sapiens.GRCh38.cdna.all.fa',
        'cand_pool_fpn': offset + 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt',
        'cdna_fp': offset + 'data/omics/genomics/fasta/cdna/GRCh38/',
    }
    simu_params = {
        'init_seq_setting': {
            'cdna_num': 50,
            'umi_unit_pattern': 3,
            'umi_unit_len': 12,
            'is_seed': True,
            'primer_len': 20,
            'cand_pool_fpn': DEFINE['cand_pool_fpn'],
            'cdna_fp': DEFINE['cdna_fp'],

            'is_sv_umi_lib': True,
            'is_sv_seq_lib': True,
            'is_sv_primer_lib': True,
            'working_dir': to('data/cdna/trimer/'),
            'umi_lib_fpn': to('data/cdna/trimer/umi.txt'),
            'seq_lib_fpn': to('data/cdna/trimer/seq.txt'),
            'primer_lib_fpn': to('data/cdna/trimer/primer.txt'),
        },
        'ampl_rate': 0.85,
        'pcr_num': 10,
        'pcr_error': 1e-3,
        'seq_error': 1e-3,
        'err_num_met': 'nbionom',
        'use_seed': False,
        'seed': None,
        'write': {
            'fastq_fp': to('data/cdna/trimer/'),
            'fastq_fn': 'simu',
        }
    }
    p = cdna(simu_params)
    print(p.ondemand())