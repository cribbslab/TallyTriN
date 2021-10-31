__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import time
import numpy as np
import pandas as pd
from simreadflow.util.random.Number import number as rannum
from simreadflow.util.sequence.symbol.Single import single as dnasgl


class error(object):

    def __init__(self, method):
        self.method = method

    def __call__(self, deal):
        from functools import wraps
        if self.method == 'default':
            func = self.postable
        @wraps(deal)
        def indexing(ph, *args, **kwargs):
            print('===>sequencing...')
            res2p = deal(ph, **kwargs)
            # print(func(res2p))
            return func(res2p)
        return indexing

    def postable(self, res2p, kind='index_by_same_len'):
        # print(res2p.keys())
        seq_stime = time.time()
        data_seq = pd.DataFrame(res2p['data_spl'], columns=['read', 'sam_id', 'source'])
        del res2p['data']
        del res2p['data_spl']
        print('======>{} reads to be sequenced'.format(data_seq.shape[0]))
        print('======>constructing the position table starts...')
        pcr_postable_stime = time.time()
        if kind == 'index_by_same_len':
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(seq_len=len(data_seq['read'][0]), num_seq=data_seq.shape[0])
        elif kind == 'index_by_lambda':
            seq_ids = []
            seq_pos_ids = []
            data_seq.apply(lambda x: self.postableLambda(x, seq_ids, seq_pos_ids), axis=1)
        else:
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(seq_len=len(data_seq['read'][0]), num_seq=data_seq.shape[0])
        pos_table = {'seq_ids': seq_ids, 'seq_pos_ids': seq_pos_ids}
        pcr_postable_etime = time.time()
        print('======>time for constructing the position table  {time:.3f}s'.format(time=pcr_postable_etime - pcr_postable_stime))
        seq_nt_num = len(seq_ids)
        print('======>{} nucleotides to be sequenced'.format(seq_nt_num))
        print('======>determining PCR error numbers starts...')
        seq_err_num_simu_stime = time.time()
        if res2p['err_num_met'] == 'binomial':
            seq_err_num = rannum().binomial(n=seq_nt_num, p=res2p['seq_error'], use_seed=True, seed=1)
        elif res2p['err_num_met'] == 'nbinomial':
            seq_err_num = rannum().nbinomial(
                n=seq_nt_num * (1 - res2p['seq_error']),
                p=1 - res2p['seq_error'],
                use_seed=True,
                seed=1
            )
        else:
            seq_err_num = rannum().binomial(n=seq_nt_num, p=res2p['seq_error'], use_seed=True, seed=1)
        print('======>{} nucleotides to be erroneous in sequencing'.format(seq_err_num))
        err_lin_ids = rannum().uniform(low=0, high=seq_nt_num, num=seq_err_num, use_seed=True, seed=1)
        # print(err_lin_ids)
        arr_err_pos = []# [[row1, col1], [row2, col2], ...]
        for i in err_lin_ids:
            arr_err_pos.append([pos_table['seq_ids'][i], pos_table['seq_pos_ids'][i]])
        pseudo_nums = rannum().uniform(low=0, high=3, num=seq_err_num, use_seed=False)
        # print(pseudo_nums)
        seq_err_num_simu_etime = time.time()
        print('======>time for determining sequencing error numbers  {time:.3f}s'.format(time=seq_err_num_simu_etime - seq_err_num_simu_stime))
        print('======>assigning sequencing errors starts...')
        seq_err_assign_stime = time.time()
        data_seq['read'] = data_seq.apply(lambda x: list(x['read']), axis=1)
        for pos_err, pseudo_num in zip(arr_err_pos, pseudo_nums):
            pcr_err_base = data_seq.loc[pos_err[0], 'read'][pos_err[1]]
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=pcr_err_base,
                    universal=True,
                ),
                reverse=True,
            )
            # print('before', data_seq.loc[pos_err[0], 'read'][pos_err[1]])
            data_seq.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # print('after', data_seq.loc[pos_err[0], 'read'][pos_err[1]])
        del arr_err_pos
        del pseudo_nums
        seq_err_assign_etime = time.time()
        print('======>time for assigning sequencing errors {time:.2f}s'.format(time=seq_err_assign_etime - seq_err_assign_stime))
        data_seq['read'] = data_seq.apply(lambda x: ''.join(x['read']), axis=1)
        res2p['data'] = data_seq.values
        del data_seq
        seq_etime = time.time()
        print('======>sequencing time: {time:.2f}s'.format(time=seq_etime - seq_stime))
        return res2p

    def postableIndexBySameLen(self, seq_len, num_seq):
        nt_ids = [i for i in range(seq_len)]
        seq_pos_ids = nt_ids * num_seq
        seq_ids = np.array([[i] * seq_len for i in range(num_seq)]).ravel().tolist()
        return seq_pos_ids, seq_ids

    def postableLambda(self, x, ids, pos_ids):
        l = len(x['read'])
        for i in range(l):
            pos_ids.append(i)
        for i in [x.name] * l:
            ids.append(i)
        return ids, pos_ids

    def epos(self, err_lin_ids, t):
        pos_err = []
        for i in err_lin_ids:
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