__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
from simreadflow.simulate.initiator.General import general as simuip
from simreadflow.pcr.Amplify import amplify as pcr
from simreadflow.sequencing.Calling import calling as seq
from simreadflow.util.sequence.fastq.Write import write as wfastq
from simreadflow.util.file.read.Reader import reader as pfreader
from simreadflow.util.random.Number import number as rannum
from simreadflow.util.sequence.symbol.Single import single as dnasgl
from Path import to


class general(object):

    def __init__(self, *args, **kwargs):
        self.args = args[0]
        self.kwargs = kwargs
        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq

    def ondemandSeqErrs(self, ):
        # ### /*** block. Init a pool of sequences ***/
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

        # ### /*** block. PCR amplification ***/
        print('->PCR amplification has started...')
        pcr_params = {
            'data': np.array(init_seqs),
            'ampl_rate': self.args['ampl_rate'],
            'err_route': self.args['err_route'],
            'pcr_error': self.args['pcr_error'],
            'pcr_num': self.args['pcr_num'],
            'err_num_met': self.args['err_num_met'],
            'use_seed': self.args['use_seed'],
            'seed': self.args['seed'],
            'recorder_nucleotide_num': [],
            'recorder_pcr_err_num': [],
            'recorder_pcr_read_num': [],
        }
        if pcr_params['err_route'] == 'minnow':
            def calclen(a):
                return len(a)
            vfunc = np.vectorize(calclen)
            pcr_params['data'] = np.hstack((vfunc(pcr_params['data'][:, 0])[:, np.newaxis], pcr_params['data'][:, 1:3]))
            print(pcr_params['data'])
            col_0 = np.array([[1] for _ in range(pcr_params['data'].shape[0])])
            cc = np.hstack((col_0, col_0))
            col_2 = pcr_params['data'][:, 1].astype(np.str)[:, np.newaxis]
            print(col_2)
            cc = np.hstack((cc, col_2))
            print(cc)

            pcr_params['mut_info'] = cc
            # pcr_params['mut_info'] = np.empty(shape=[0, 3])

        pcr = self.pcr(pcr_params=pcr_params).np()
        print(pcr.keys())
        print('->PCR amplification has completed.')

        if pcr_params['err_route'] == 'minnow':
            pcr['data'] = self.subsamplingMinnow(pcr_dict=pcr)

        # ### /*** block. Sequencing ***/
        # print('->Sequencing has started...')
        # for id, iseq_err in enumerate(self.args['seq_errors']):
        #     seq_params = {
        #         'data': pcr['data'],
        #         'seq_sub_spl_rate': self.args['seq_sub_spl_rate'],
        #         'seq_error': iseq_err,
        #         'err_num_met': self.args['err_num_met'],
        #         'use_seed': self.args['use_seed'],
        #         'seed': self.args['seed'],
        #     }
        #     seq = self.seq(seq_params=seq_params).np()
        #     print('->Sequencing has completed.')
        #     print('->Write seqs in fastq format')
        #     self.wfastq().togz(
        #         list_2d=seq['data'],
        #         sv_fp=self.args['write']['fastq_fp'],
        #         fn=self.args['write']['fastq_fn'] + 'seq_err_' + str(id),
        #         symbol='-',
        #     )
        #     del seq
        return

    def subsamplingMinnow(self, pcr_dict):
        umi_map = pfreader().generic(self.args['init_seq_setting']['umi_lib_fpn'])[0].to_dict()
        # print(umi_map)
        nn = pcr_dict['data'].shape[0]
        spl_ids = rannum().uniform(
            low=0, high=nn, num=500, use_seed=False, seed=1
        )
        # print(spl_ids)
        spl_id_map = self.tactic6(pcr_dict['data'][:, [1, 2]])
        spl_mut_info = pcr_dict['mut_info'][spl_ids]
        keys = spl_mut_info[:, 2]
        print(keys)
        pos_dict = self.tactic6(pcr_dict['mut_info'][:, [2, 0]])
        base_dict = self.tactic6(pcr_dict['mut_info'][:, [2, 1]])
        # print()
        # print(self.tactic6(base_np))
        res_data = []
        for key in keys:
            mol_id = key.split('_')[0]
            k = key.split('_')[1:]
            print('kkk', key, k)
            read = umi_map[int(mol_id)]
            for i in range(len(k)):
                print('id', i)
                sub_k = mol_id + '_' + '_'.join(k[:i+1]) if k != [] else mol_id
                print(pos_dict[sub_k], base_dict[sub_k])
                read = self.change(read, pos_list=pos_dict[sub_k], base_list=base_dict[sub_k])
                print(read)
            print(read)
            res_data.append([
                read,  # read
                str(mol_id) + '_' + '_'.join(k) if k != [] else str(mol_id),  # sam id
                spl_id_map[str(mol_id) + '_' + '_'.join(k)] if k != [] else 'init',  # source
            ])
            print(read)
        print(np.array(res_data).shape)
        return np.array(res_data)

    def change(self, read, pos_list, base_list):
        read_l = list(read)
        for i, pos in enumerate(pos_list):
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=read_l[pos],
                    universal=True,
                ),
                reverse=True,
            )
            read_l[pos] = dna_map[base_list[i]]
        return ''.join(read_l)

    def findListDuplicates(self, l):
        seen = set()
        seen_add = seen.add
        seen_twice = set(x for x in l if x in seen or seen_add(x))
        return list(seen_twice)

    def tactic6(self, arr_2d):
        result = {}
        len_arr = len(arr_2d[0])
        if len_arr == 2:
            for item in arr_2d:
                result[item[0]] = item[1]
        else:
            for item in arr_2d:
                result[item[0]] = item[1:]
        return result

    def ondemandPCRErrs(self, ):
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
        for id, ipcr_err in enumerate(self.args['pcr_errors']):
            pcr_params = {
                'data': np.array(init_seqs),
                'ampl_rate': self.args['ampl_rate'],
                'pcr_error': ipcr_err,
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
            seq_params = {
                'data': pcr['data'],
                'seq_sub_spl_rate': self.args['seq_sub_spl_rate'],
                'seq_error': self.args['seq_error'],
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
                fn=self.args['write']['fastq_fn'] + 'pcr_err_' + str(id),
                symbol='-',
            )
            del seq
        return

    def ondemandUMILens(self, ):
        # /*** block. Init a pool of sequences ***/
        # print(self.args)
        # print(self.args['init_seq_setting']['umi_unit_lens'])
        for id, iumi_len in enumerate(self.args['init_seq_setting']['umi_unit_lens']):
            print('->Init pool of sequences has started.')
            init_seqs = simuip(
                seq_num=self.args['init_seq_setting']['seq_num'],
                umi_unit_pattern=self.args['init_seq_setting']['umi_unit_pattern'],
                umi_unit_len=iumi_len,
                is_seed=self.args['init_seq_setting']['is_seed'],
                is_sv_umi_lib=self.args['init_seq_setting']['is_sv_umi_lib'],
                umi_lib_fpn=self.args['init_seq_setting']['umi_lib_fp'] + '/umi_' + str(iumi_len) + '.txt',
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
            seq_params = {
                'data': pcr['data'],
                'seq_sub_spl_rate': self.args['seq_sub_spl_rate'],
                'seq_error': self.args['seq_error'],
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
                fn=self.args['write']['fastq_fn'] + 'umi_len_' + str(iumi_len),
                symbol='-',
            )
            del seq
        return

    def ondemandAmplRates(self, ):
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
        for id, ampl_rate in enumerate(self.args['ampl_rates']):
            pcr_params = {
                'data': np.array(init_seqs),
                'ampl_rate': ampl_rate,
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
            seq_params = {
                'data': pcr['data'],
                'seq_sub_spl_rate': self.args['seq_sub_spl_rate'],
                'seq_error': self.args['seq_error'],
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
                fn=self.args['write']['fastq_fn'] + 'ampl_rate_' + str(id),
                symbol='-',
            )
            del seq
        return

    def ondemandPCRNums(self, ):
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
        for id, ipcr_num in enumerate(self.args['pcr_nums']):
            pcr_params = {
                'data': np.array(init_seqs),
                'ampl_rate': self.args['ampl_rate'],
                'pcr_error': self.args['pcr_error'],
                'pcr_num': ipcr_num,
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
            seq_params = {
                'data': pcr['data'],
                'seq_sub_spl_rate': self.args['seq_sub_spl_rate'],
                'seq_error': self.args['seq_error'],
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
                fn=self.args['write']['fastq_fn'] + 'pcr_num_' + str(ipcr_num),
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
    umi_unit_len = 6
    simu_params = {
        'init_seq_setting': {
            'seq_num': 50,
            'umi_unit_pattern': 1,
            'umi_unit_len': umi_unit_len,
            'is_seed': True,
            'is_sv_umi_lib': True,
            'working_dir': to('data/simu/monomer/umi_lens/'),
            'umi_lib_fp': to('data/simu/monomer/umi_lens/'),
            'umi_lib_fpn': to('data/simu/monomer/umi_lens/umi.txt'),
            'condis': ['umi'],
            'sim_thres': 3,
            'permutation': 1,
        },
        'ampl_rate': 0.85,
        'pcr_num': 16, # 60,000,000
        'err_num_met': 'nbinomial',
        'pcr_error': 1e-1,
        'seq_error': 1e-1,
        "ampl_rates": np.linspace(0.1, 1, 10),
        "umi_unit_lens": np.arange(6, 36 + 1, 1),
        "umi_nums": np.arange(20, 140 + 20, 20),
        "pcr_nums": np.arange(1, 20 + 1, 1),
        'pcr_errors': [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3],
        'seq_errors': [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3],
        'seq_sub_spl_rate': 0.3333,
        'use_seed': False,
        'seed': None,
        'write': {
            'fastq_fp': to('data/simu/monomer/umi_lens/'),
            'fastq_fn': '',
        }
    }
    p = general(simu_params)
    print(p.ondemandSeqErrs())
    # print(p.ondemandUMILens())
    # print(p.ondemandPCRErrs())
    # print(p.ondemandPCRNums())
    # print(p.ondemandAmplRates())