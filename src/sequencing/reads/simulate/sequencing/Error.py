__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__author__ = "Adam Cribbs lab"

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
            print('indexing...')
            res2p = deal(ph, **kwargs)
            # print(func(res2p))
            return func(res2p)
        return indexing

    def tomap(self, res2p):
        data = pd.DataFrame(res2p['data'], columns=['read', 'sam_id', 'source'])
        df_read_len = data['read'].apply(lambda x: len(x))
        data['read'] = data.apply(lambda x: list(x['read']), axis=1)
        print(data.shape)
        # print(data)
        num_bases = df_read_len.sum()
        # print(num_bases)
        seq_err_num = rannum().binomial(n=num_bases, p=res2p['seq_error'], use_seed=True, seed=1)
        # print(seq_err_num)
        err_lin_ids = rannum().uniform(low=0, high=num_bases, num=seq_err_num, use_seed=True, seed=1)
        # print(err_lin_ids)
        err_arr2d_pos = self.epos(err_lin_ids=err_lin_ids, t=df_read_len.tolist())
        print('err len: ', len(err_arr2d_pos))
        pseudo_nums = rannum().uniform(low=0, high=3, num=seq_err_num, use_seed=True, seed=1)
        # print(pseudo_nums)
        for pos_err, pseudo_num in zip(err_arr2d_pos, pseudo_nums):
            pcr_err_base = data.loc[pos_err[0], 'read'][pos_err[1]]
            dna_map = dnasgl().todict(
                bases=dnasgl().getEleTrimmed(
                    ele_loo=pcr_err_base,
                    universal=True,
                ),
                reverse=True,
            )
            # print('before', data.loc[pos_err[0], 'read'][pos_err[1]])
            data.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # print('after', data.loc[pos_err[0], 'read'][pos_err[1]])
        data['read'] = data.apply(lambda x: ''.join(x['read']), axis=1)
        res2p['data'] = data.values
        print(res2p['data'].shape)
        return res2p

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