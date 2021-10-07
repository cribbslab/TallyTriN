__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from src.sequencing.reads.simulate.initiator.GeneralCondi import generalCondi as simuip
from src.sequencing.reads.simulate.pcr.Amplify import amplify as pcr
from src.sequencing.reads.simulate.sequencing.Calling import calling as seq
from src.util.sequence.fastq.Write import write as wfastq
from Path import to


class generalCondi(object):

    def __init__(self, *args, **kwargs):
        self.args = args[0]
        self.kwargs = kwargs
        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq

    def ondemand(self, ):
        # /*** block. Init a pool of sequences ***/
        print('->Init pool of sequences has started.')
        init_seqs = simuip(
            seq_num=self.args['init_seq_setting']['seq_num'],
            seq_len=self.args['init_seq_setting']['seq_len'],
            spacer_len=self.args['init_seq_setting']['spacer_len'],
            umi_unit_pattern=self.args['init_seq_setting']['umi_unit_pattern'],
            umi_unit_len=self.args['init_seq_setting']['umi_unit_len'],
            is_seed=self.args['init_seq_setting']['is_seed'],
            is_sv_umi_lib=self.args['init_seq_setting']['is_sv_umi_lib'],
            is_sv_spacer_lib=self.args['init_seq_setting']['is_sv_spacer_lib'],
            is_sv_seq_lib=self.args['init_seq_setting']['is_sv_seq_lib'],
            umi_lib_fpn=self.args['init_seq_setting']['umi_lib_fpn'],
            spacer_lib_fpn=self.args['init_seq_setting']['spacer_lib_fpn'],
            seq_lib_fpn=self.args['init_seq_setting']['seq_lib_fpn'],
            working_dir=self.args['init_seq_setting']['working_dir'],
            condis=self.args['init_seq_setting']['condis'],
            permutation=self.args['init_seq_setting']['permutation'],
        ).generate()
        self.args['init_seqs'] = init_seqs
        print('->Init pool of sequences has completed.')
        # print(init_seqs)
        # /*** block. PCR amplification ***/
        print('->PCR amplification has started...')
        pcr = self.pcr(self.args).np()
        print('->PCR amplification has completed.')
        # /*** block. Sequencing ***/
        print(pcr.keys())
        print('->Sequencing has started...')
        pcr['seq_error'] = self.args['seq_error']
        seq = self.seq(pcr).np()
        print('->Sequencing has completed.')
        print('->Write seqs in fastq format')
        self.wfastq().togz(
            list_2d=seq['data'],
            sv_fp=self.args['write']['fastq_fp'],
            fn=self.args['write']['fastq_fn'],
            symbol='-',
        )
        return


if __name__ == "__main__":
    dis = '../../../'
    offset = '../' * 4 + dis
    DEFINE = {
        '': '',
    }
    umi_unit_len = 12
    simu_params = {
        'init_seq_setting': {
            'seq_num': 100,
            'umi_unit_pattern': 1,
            'umi_unit_len': umi_unit_len,
            'is_seed': True,
            'seq_len': 100 - umi_unit_len,
            'is_sv_umi_lib': True,
            'is_sv_seq_lib': True,
            'working_dir': to('data/simu/'),
            'umi_lib_fpn': to('data/simu/umi.txt'),
            'seq_lib_fpn': to('data/simu/seq.txt'),
        },
        'ampl_rate': 0.85,
        'pcr_num': 3,
        'err_num_met': 'bionom',
        'pcr_error': 1e-3,
        'seq_error': 1e-3,
        'use_seed': True,
        'seed': None,
        'write': {
            'fastq_fp': to('data/'),
            'fastq_fn': 'simu',
        }
    }
    p = generalCondi(simu_params)
    print(p.ondemand())