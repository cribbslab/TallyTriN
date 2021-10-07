__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
from src.sequencing.reads.simulate.initiator.UMI import umi as simuip
from src.sequencing.reads.simulate.pcr.Amplify import amplify as pcr
from src.sequencing.reads.simulate.sequencing.Calling import calling as seq
from src.util.sequence.fastq.Write import write as wfastq
from Path import to


class umi(object):

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
            umi_unit_pattern=self.args['init_seq_setting']['umi_unit_pattern'],
            umi_unit_len=self.args['init_seq_setting']['umi_unit_len'],
            is_seed=self.args['init_seq_setting']['is_seed'],
            is_sv_umi_lib=self.args['init_seq_setting']['is_sv_umi_lib'],
            umi_lib_fpn=self.args['init_seq_setting']['umi_lib_fpn'],
            working_dir=self.args['init_seq_setting']['working_dir'],
            condis=self.args['init_seq_setting']['condis'],
            permutation=self.args['init_seq_setting']['permutation'],
        ).pooling()
        print('->Init pool of sequences has completed.')
        # print(init_seqs)

        # /*** block. PCR amplification ***/
        print('->PCR amplification has started...')

        pcr_params = {
            'data': np.array(init_seqs),
            'ampl_rate': self.args['ampl_rate'],
            'pcr_error': self.args['pcr_error'],
            'pcr_num': self.args['pcr_num'],
            'err_num_met': self.args['err_num_met'],
            'use_seed': self.args['use_seed'],
            'seed': self.args['seed'],
            'recorder_nucleotide_num': [],
            'recorder_pcr_err_num': [],
            'recorder_pcr_read_num': [],
        }
        pcr = self.pcr(pcr_params=pcr_params).np()
        print(pcr.keys())
        print('->PCR amplification has completed.')

        # /*** block. Sequencing ***/
        print('->Sequencing has started...')
        for id, iseq_err in enumerate(self.args['seq_error']):
            seq_params = {
                'data': pcr['data'],
                'seq_sub_spl_rate': self.args['seq_sub_spl_rate'],
                'seq_error': iseq_err,
                'err_num_met': self.args['err_num_met'],
                'use_seed': self.args['use_seed'],
                'seed': self.args['seed'],
            }
            seq = self.seq(seq_params=seq_params).np()
            print('->Sequencing has completed.')
            print('->Write seqs in fastq format')
            self.wfastq().togz(
                list_2d=seq['data'],
                sv_fp=self.args['write']['fastq_fp'],
                fn=self.args['write']['fastq_fn'] + 'seq_err_' + str(id),
                symbol='-',
            )
            del seq
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
            'seq_num': 1000,
            'umi_unit_pattern': 1,
            'umi_unit_len': umi_unit_len,
            'is_seed': True,
            'is_sv_umi_lib': True,
            'working_dir': to('data/simu/umi/seq_errs/monomer/'),
            'umi_lib_fpn': to('data/simu/umi/seq_errs/monomer/umi.txt'),
            'condis': ['umi'],
            'sim_thres': 3,
            'permutation': 1,
        },
        'ampl_rate': 0.85,
        'pcr_num': 16, # 60,000,000
        'err_num_met': 'nbinomial',
        'pcr_error': 1e-4,
        # 'seq_error': 1e-1,
        'seq_error': [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3],
        'seq_sub_spl_rate': 0.3333,
        'use_seed': False,
        'seed': None,
        'write': {
            'fastq_fp': to('data/simu/umi/seq_errs/monomer/'),
            'fastq_fn': '',
        }
    }
    p = umi(simu_params)
    print(p.ondemand())