__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__author__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd
from src.util.random.Number import number as rannum


class error(object):

    def __init__(self, method):
        self.method = method

    def __call__(self, deal):
        from functools import wraps
        if self.method == 'tomap':
            func = self.tomap
        @wraps(deal)
        def indexing(ph, *args, **kwargs):
            print('indexing...')
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
            # t = np.arange(len(d))[:, np.newaxis]
            # map_tab = np.concatenate((t, d), axis=1)
            # print(self.tactic1(map_tab))

        :param res2p:
        :return:
        """
        data_pcr = pd.DataFrame(res2p['data_spl'], columns=['read', 'sam_id', 'source'])
        res2p['recorder_num_pcr_reads'].append(data_pcr.shape[0])
        ampl_read_len_df = data_pcr['read'].apply(lambda x: len(x))
        data_pcr['read'] = data_pcr.apply(lambda x: list(x['read']), axis=1)
        # print(data_pcr.shape)
        # print(data_pcr)
        num_bases = ampl_read_len_df.sum()
        # print(num_bases)
        res2p['recorder_num_bases'].append(num_bases)
        pcr_err_num = rannum().binomial(n=num_bases, p=res2p['pcr_error'], use_seed=True, seed=1)
        # print(pcr_err_num)
        res2p['recorder_num_pcr_errs'].append(pcr_err_num)
        spl_ids = rannum().uniform(low=0, high=num_bases, num=pcr_err_num, use_seed=True, seed=1)
        # print(spl_ids)
        arr_err_pos = self.epos(spl_ids=spl_ids, t=ampl_read_len_df.tolist())
        # print(arr_err_pos)
        pseudo_nums = rannum().uniform(low=0, high=3, num=pcr_err_num, use_seed=True, seed=1)
        # print(pseudo_nums)
        from src.util.sequence.symbol.Single import single as dnasgl
        for pos_err, pseudo_num in zip(arr_err_pos, pseudo_nums):
            pcr_err_base = data_pcr.loc[pos_err[0], 'read'][pos_err[1]]
            dna_map = dnasgl().todict(
                bases=dnasgl().getEleTrimmed(
                    ele_loo=pcr_err_base,
                    universal=True,
                ),
                reverse=True,
            )
            # print('before', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
            data_pcr.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # print('after', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
        data_pcr['read'] = data_pcr.apply(lambda x: ''.join(x['read']), axis=1)
        data_pcr['source'] = 'pcr_' + str(res2p['ipcr'] + 1)
        # print(data_pcr.values)
        data_pcr = np.array(data_pcr)
        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        print(res2p['data'].shape)
        print(res2p['recorder_num_pcr_reads'])
        print(res2p['recorder_num_bases'])
        print(res2p['recorder_num_pcr_errs'])
        return res2p

    def epos(self, spl_ids, t):
        pos_err = []
        for i in spl_ids:
            cc = 0
            for id, j in enumerate(t):
                cc += j
                if cc >= i:
                    # print(cc, i)
                    pos_err.append([id, (j - 1) - (cc % i)])
                    break
        return pos_err

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

    def todict(self, bases, reverse=False):
        aa_dict = {}
        for k, v in enumerate(bases):
            aa_dict[v] = k
        if reverse:
            aa_dict = {v: k for k, v in aa_dict.items()}
        return aa_dict