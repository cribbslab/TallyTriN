__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__author__ = "Adam Cribbs lab"

import time
import numpy as np
import pandas as pd
from src.util.random.Number import number as rannum
from src.util.sequence.symbol.Single import single as dnasgl


class error(object):

    def __init__(self, method):
        self.method = method

    def __call__(self, deal):
        from functools import wraps
        if self.method == 'tomap':
            func = self.tomap
        @wraps(deal)
        def indexing(ph, *args, **kwargs):
            print('------>error making...')
            res2p = deal(ph, **kwargs)
            # print(func(res2p))
            return func(res2p)
        return indexing

    def tomap(self, res2p):
        """
        ..  @test:
            ------
            # dna = ['A', 'T', 'G', 'C']
            # dna.remove(pcr_err_base)
            # dna_map = self.todict(dna, reverse=True)
            # data_pcr.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # pcr_err_base = data_pcr[pos_err[0]][0][pos_err[1]]

            # d = [
            #     [i, j] for i in ampl_read_len_df.index for j in np.arange(ampl_read_len_df[i])
            # ]
            # print(len(d))
            # ampl_read_len_list = np.arange(len(d))[:, np.newaxis]
            # map_tab = np.concatenate((ampl_read_len_list, d), axis=1)
            # print(self.tactic1(map_tab))

        :param res2p:
        :return:
        """
        # print(res2p.keys())
        pcr_stime = time.time()
        data_pcr = pd.DataFrame(res2p['data_spl'], columns=['read', 'sam_id', 'source'])
        res2p['recorder_pcr_read_num'].append(data_pcr.shape[0])
        print('------>{} reads to be amplified'.format(data_pcr.shape[0]))
        print('------>constructing the position table starts...')
        pcr_postable_stime = time.time()
        seq_ids = []
        seq_pos_ids = []
        data_pcr.apply(lambda x: self.postable(x, seq_ids, seq_pos_ids), axis=1)
        pos_table = {'seq_ids': seq_ids, 'seq_pos_ids': seq_pos_ids}
        pcr_postable_etime = time.time()
        print('------>time for constructing the position table  {time:.3f}s'.format(time=pcr_postable_etime - pcr_postable_stime))
        ampl_nt_num = len(seq_ids)
        data_pcr['read'] = data_pcr.apply(lambda x: list(x['read']), axis=1)
        print('------>{} nucleotides to be amplified'.format(ampl_nt_num))
        res2p['recorder_nucleotide_num'].append(ampl_nt_num)
        print('------>determining PCR error numbers starts...')
        pcr_err_num_simu_stime = time.time()
        if res2p['err_num_met'] == 'bionom':
            pcr_err_num = rannum().binomial(n=ampl_nt_num, p=res2p['pcr_error'], use_seed=True, seed=res2p['ipcr'] + 1)
        elif res2p['err_num_met'] == 'nbionom':
            pcr_err_num = rannum().nbinomial(
                n=ampl_nt_num*(1-res2p['pcr_error']),
                p=1-res2p['pcr_error'],
                use_seed=True,
                seed=res2p['ipcr'] + 1
            )
        else:
            pcr_err_num = rannum().binomial(n=ampl_nt_num, p=res2p['pcr_error'], use_seed=True, seed=res2p['ipcr'] + 1)
        print('------>{} nucleotides to be erroneous at this PCR'.format(pcr_err_num))
        res2p['recorder_pcr_err_num'].append(pcr_err_num)
        spl_nt_ids = rannum().uniform(low=0, high=ampl_nt_num, num=pcr_err_num, use_seed=True, seed=res2p['ipcr'] + 1)
        arr_err_pos = []
        for i in spl_nt_ids:
            arr_err_pos.append([pos_table['seq_ids'][i], pos_table['seq_pos_ids'][i]])
        pseudo_nums = rannum().uniform(low=0, high=3, num=pcr_err_num, use_seed=False)
        # print(pseudo_nums)
        pcr_err_num_simu_etime = time.time()
        print('------>time for determining PCR error numbers  {time:.3f}s'.format(time=pcr_err_num_simu_etime - pcr_err_num_simu_stime))
        print('------>assigning PCR errors starts...')
        pcr_err_assign_stime = time.time()
        for pos_err, pseudo_num in zip(arr_err_pos, pseudo_nums):
            pcr_err_base = data_pcr.loc[pos_err[0], 'read'][pos_err[1]]
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=pcr_err_base,
                    universal=True,
                ),
                reverse=True,
            )
            # print('before', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
            data_pcr.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # print('after', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
        pcr_err_assign_etime = time.time()
        print('------>time for assigning PCR errors {time:.2f}s'.format(time=pcr_err_assign_etime - pcr_err_assign_stime))
        print('------>merging the PCR duplicates and all previous sequences starts...')
        pcr_merge_stime = time.time()
        data_pcr['read'] = data_pcr.apply(lambda x: ''.join(x['read']), axis=1)
        data_pcr['source'] = 'pcr-' + str(res2p['ipcr'] + 1)
        # print(data_pcr.values)
        data_pcr = np.array(data_pcr)
        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        pcr_merge_etime = time.time()
        print('------>time for merging sequences {time:.2f}s'.format(time=pcr_merge_etime - pcr_merge_stime))
        pcr_etime = time.time()
        print('------>Summary report:')
        print('--------->PCR time: {}'.format(pcr_etime - pcr_stime))
        print('--------->the dimensions of the data: number of reads: {}'.format(res2p['data'].shape))
        print('--------->the number of reads at this PCR: {}, '.format(res2p['recorder_pcr_read_num']))
        print('--------->the number of nucleotides at this PCR: {}, '.format(res2p['recorder_nucleotide_num']))
        print('--------->the number of errors at this PCR: {}, '.format(res2p['recorder_pcr_err_num']))
        return res2p

    def postable(self, x, ids, pos_ids):
        l = len(x['read'])
        for i in range(l):
            pos_ids.append(i)
        for i in [x.name] * l:
            ids.append(i)
        return ids, pos_ids

    def epos_deprecate(self, spl_nt_ids, ampl_read_len_list):
        """
        ..  @call:
            ------
            arr_err_pos = self.epos(spl_nt_ids=spl_nt_ids, ampl_read_len_list=ampl_read_len_list)

        :param spl_nt_ids: is a list of ids of nucleotides to be modified.
        :param ampl_read_len_list: a list of lengths of reads to be amplified.
        :return:
        """
        pos_err = []
        for i in spl_nt_ids:
            cc = 0
            for id, j in enumerate(ampl_read_len_list):
                cc += j
                if cc >= i:
                    # print(cc, i)
                    pos_err.append([id, (j - 1) - (cc % i)])
                    break
        return pos_err

    def epos2_deprecate(self, ampl_read_len_list):
        ids = []
        pos = []
        for k, e in enumerate(ampl_read_len_list):
            l = [k] * e
            ids = ids + l
            p = [i for i in range(e)]
            pos = pos + p
        return np.array([pos, ids])

    def tactic1(self, arr_2d):
        result = {}
        len_arr = len(arr_2d[0])
        if len_arr == 2:
            for item in arr_2d:
                result[item[0]] = item[1]
        else:
            for item in arr_2d:
                result[item[0]] = item[1:]
        return result

    def todict(self, nucleotides, reverse=False):
        aa_dict = {}
        for k, v in enumerate(nucleotides):
            aa_dict[v] = k
        if reverse:
            aa_dict = {v: k for k, v in aa_dict.items()}
        return aa_dict