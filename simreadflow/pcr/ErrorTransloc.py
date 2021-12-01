__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import time
import numpy as np
import pandas as pd
from simreadflow.util.random.Number import number as rannum
from simreadflow.util.sequence.symbol.Single import single as dnasgl


class errorTransloc(object):

    def __init__(self, method):
        self.method = method

    def __call__(self, deal):
        from functools import wraps
        if self.method == 'default':
            func = self.postable
        @wraps(deal)
        def indexing(ph, *args, **kwargs):
            print('======>error making...')
            res2p = deal(ph, **kwargs)
            # print(func(res2p))
            return func(res2p)
        return indexing

    def postable(self, res2p, kind='index_by_same_len'):
        """
        ..  test:
            -----
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

        ..  Method:
            -------
            =>string concat
                for pos_err, pseudo_num in zip(arr_err_pos, pseudo_nums):
                    tar_read = data_pcr.loc[pos_err[0], 'read']
                    pcr_err_base = tar_read[pos_err[1]]
                    dna_map = dnasgl().todict(
                        nucleotides=dnasgl().getEleTrimmed(
                            ele_loo=pcr_err_base,
                            universal=True,
                        ),
                        reverse=True,
                    )
                    # print('before', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
                    data_pcr.loc[pos_err[0], 'read'] = tar_read[:pos_err[1]] + dna_map[pseudo_num] + tar_read[pos_err[1]+1:]
        :param res2p:
        :return:
        """
        # print(res2p.keys())
        pcr_stime = time.time()
        data_pcr = pd.DataFrame(res2p['data_spl'], columns=['read', 'r1_id', 'r2_id', 'transloc_stat', 'transloc_side', 'sam_id', 'source'])
        data_pcr['transloc_stat'] = 'fake_no'
        transloc_num = rannum().binomial(n=data_pcr.shape[0], p=res2p['transloc_rate'], use_seed=True, seed=res2p['ipcr'] + 1)
        # print('++++++++> {}'.format(transloc_num))
        transloc_rid = rannum().choice(high=data_pcr.shape[0], num=transloc_num*2, replace=False)
        transloc_rids = np.reshape(transloc_rid, (transloc_num, 2))
        # print('++++++++> {}'.format(transloc_rid))
        # print('++++++++> {}'.format(transloc_rids.shape))

        for pos_read in transloc_rids:
            transloc_read_before_1st = data_pcr.loc[pos_read[0]]
            transloc_read_before_2nd = data_pcr.loc[pos_read[1]]
            # print('1st before {}'.format(transloc_read_before_1st['read']))
            # print(transloc_read_before_1st['r1_id'])
            # print(transloc_read_before_1st['r2_id'])
            # print(transloc_read_before_1st['transloc_stat'])
            # print(transloc_read_before_1st['transloc_side'])
            # print('2nd before {}'.format(transloc_read_before_2nd['read']))
            # print(transloc_read_before_2nd['r1_id'])
            # print(transloc_read_before_2nd['r2_id'])
            # print(transloc_read_before_2nd['transloc_stat'])
            # print(transloc_read_before_2nd['transloc_side'])
            transloc_read_after_1st = transloc_read_before_1st['read'][:36] + transloc_read_before_2nd['read'][36:]
            transloc_read_after_2nd = transloc_read_before_2nd['read'][:36] + transloc_read_before_1st['read'][36:]
            # print('1st after {}'.format(transloc_read_after_1st))
            # print('2nd after {}'.format(transloc_read_after_2nd))
            data_pcr.loc[pos_read[0], 'read'] = transloc_read_after_1st
            data_pcr.loc[pos_read[1], 'read'] = transloc_read_after_2nd
            data_pcr.loc[pos_read[0], 'transloc_stat'] = 'fake_yes'
            data_pcr.loc[pos_read[1], 'transloc_stat'] = 'fake_yes'
            data_pcr.loc[pos_read[0], 'transloc_side'] = 'r'
            data_pcr.loc[pos_read[1], 'transloc_side'] = 'r'
            tmp_r2_id_1st = data_pcr.loc[pos_read[0], 'r2_id']
            tmp_r2_id_2nd = data_pcr.loc[pos_read[1], 'r2_id']
            data_pcr.loc[pos_read[0], 'r2_id'] = tmp_r2_id_2nd
            data_pcr.loc[pos_read[1], 'r2_id'] = tmp_r2_id_1st
            # print(transloc_read_before_1st['transloc_stat'])
            # print(transloc_read_before_2nd['transloc_stat'])
            # print(data_pcr.loc[pos_read[0], 'read'])
            # print(data_pcr.loc[pos_read[1], 'read'])
            # print(data_pcr.loc[pos_read[0], 'r2_id'])
            # print(data_pcr.loc[pos_read[1], 'r2_id'])
        print(data_pcr['transloc_stat'])
        del res2p['data_spl']
        # print(len(data_pcr['read'][0]))
        # print(data_pcr.shape[0])
        res2p['recorder_pcr_read_num'].append(data_pcr.shape[0])
        print('======>{} reads to be amplified'.format(data_pcr.shape[0]))
        print('======>constructing the position table starts...')
        # print(data_pcr)
        pcr_postable_stime = time.time()
        if kind == 'index_by_same_len':
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(seq_len=len(data_pcr['read'][0]), num_seq=data_pcr.shape[0])
        elif kind == 'index_by_lambda':
            seq_ids = []
            seq_pos_ids = []
            data_pcr.apply(lambda x: self.postableLambda(x, seq_ids, seq_pos_ids), axis=1)
        else:
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(seq_len=len(data_pcr['read'][0]), num_seq=data_pcr.shape[0])
        # print(seq_ids)
        pos_table = {'seq_ids': seq_ids, 'seq_pos_ids': seq_pos_ids}
        print('======>time for constructing the position table  {time:.3f}s'.format(time=time.time() - pcr_postable_stime))
        ampl_nt_num = len(seq_ids)
        print('======>{} nucleotides to be amplified'.format(ampl_nt_num))
        res2p['recorder_nucleotide_num'].append(ampl_nt_num)
        print('======>determining PCR error numbers starts...')
        pcr_err_num_simu_stime = time.time()
        if res2p['err_num_met'] == 'binomial':
            pcr_err_num = rannum().binomial(n=ampl_nt_num, p=res2p['pcr_error'], use_seed=True, seed=res2p['ipcr'] + 1)
        elif res2p['err_num_met'] == 'nbinomial':
            pcr_err_num = rannum().nbinomial(
                n=ampl_nt_num*(1-res2p['pcr_error']),
                p=1-res2p['pcr_error'],
                use_seed=True,
                seed=res2p['ipcr'] + 1
            )
            # print('++++++++> {}'.format(pcr_err_num))
        else:
            pcr_err_num = rannum().binomial(n=ampl_nt_num, p=res2p['pcr_error'], use_seed=True, seed=res2p['ipcr'] + 1)
        print('======>{} nucleotides to be erroneous at this PCR'.format(pcr_err_num))
        res2p['recorder_pcr_err_num'].append(pcr_err_num)
        spl_nt_ids = rannum().uniform(low=0, high=ampl_nt_num, num=pcr_err_num, use_seed=True, seed=res2p['ipcr'] + 1)
        arr_err_pos = []
        for i in spl_nt_ids:
            arr_err_pos.append([pos_table['seq_ids'][i], pos_table['seq_pos_ids'][i]])
        pseudo_nums = rannum().uniform(low=0, high=3, num=pcr_err_num, use_seed=False)
        # print(pseudo_nums)
        print('======>time for determining PCR error numbers  {time:.3f}s'.format(time=time.time() - pcr_err_num_simu_stime))
        print('======>assigning PCR errors starts...')
        pcr_err_assign_stime = time.time()
        data_pcr['read'] = data_pcr.apply(lambda x: list(x['read']), axis=1)
        # caser = self.todict5(arr_err_pos)
        # data_pcr['read'] = data_pcr.apply(lambda x: self.worktable(x, caser), axis=1)
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
        print('======>time for assigning PCR errors {time:.2f}s'.format(time=time.time() - pcr_err_assign_stime))
        print('======>merging the PCR duplicates and all previous sequences starts...')
        del arr_err_pos
        del spl_nt_ids
        pcr_merge_stime = time.time()
        data_pcr['read'] = data_pcr.apply(lambda x: ''.join(x['read']), axis=1)
        data_pcr['source'] = 'pcr' + str(res2p['ipcr'] + 1)
        # print(data_pcr.values)
        data_pcr = np.array(data_pcr)
        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        del data_pcr
        print('======>time for merging sequences {time:.2f}s'.format(time=time.time() - pcr_merge_stime))
        print('======>Summary report:')
        print('=========>PCR time: {time:.2f}s'.format(time=time.time() - pcr_stime))
        print('=========>the dimensions of the data: number of reads: {}'.format(res2p['data'].shape))
        print('=========>the number of reads at this PCR: {}, '.format(res2p['recorder_pcr_read_num']))
        print('=========>the number of nucleotides at this PCR: {}, '.format(res2p['recorder_nucleotide_num']))
        print('=========>the number of errors at this PCR: {}, '.format(res2p['recorder_pcr_err_num']))
        return res2p

    def postableIndexBySameLen(self, seq_len, num_seq):
        nt_ids = [i for i in range(seq_len)]
        seq_pos_ids = nt_ids * num_seq
        seq_ids = np.array([[i] * seq_len for i in range(num_seq)]).ravel().tolist()
        return seq_pos_ids, seq_ids

    def postableLambda(self, x, ids, pos_ids):
        for i in range(x['read_len']):
            pos_ids.append(i)
        ids = ids + [x.name] * x['read_len']
        return ids, pos_ids

    def worktable(self, x, ids):
        if x.name not in ids:
            return x['read']
        else:
            tt = x['read']
            for i in ids[x.name]:
                pcr_err_base = tt[i]
                dna_map = dnasgl().todict(
                    nucleotides=dnasgl().getEleTrimmed(
                        ele_loo=pcr_err_base,
                        universal=True,
                    ),
                    reverse=True,
                )
                tt = tt[:i] + dna_map[2] + tt[i+1:]
            return tt

    def todict(self, nucleotides, reverse=False):
        aa_dict = {}
        for k, v in enumerate(nucleotides):
            aa_dict[v] = k
        if reverse:
            aa_dict = {v: k for k, v in aa_dict.items()}
        return aa_dict

    def todict5(self, arr_2d):
        result = {}
        for item in arr_2d:
            result[item[0]] = []
        for item in arr_2d:
            if item[0] in result.keys():
                result[item[0]].append(item[1])
                # print(result[item[0]])
        return result