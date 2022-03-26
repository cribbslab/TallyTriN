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
        if pcr_params['err_route'] == 'tree':
            pcr_params['data'] = pcr_params['data'][:, 1:3]
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

        # ### /*** block. Subsampling if tree ***/
        if pcr_params['err_route'] == 'tree':
            pcr['data'] = self.subsamplingPCRTree(pcr_dict=pcr)

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

    def subsamplingPCRTree(self, pcr_dict):
        """

        Notes
        -----
            bool_flags_ = [[] for _ in range(len(uniq_mol_map_new))]
            realvalued_flags_ = [[] for _ in range(len(uniq_mol_map_new))]
            uniq_mol_map_new_ = [[] for _ in range(len(uniq_mol_map_new))]
            for ii, u in enumerate(uniq_mol_map_new):
                for jj, v in enumerate(u):
                    if v != '0':
                        bool_flags_[ii].append(bool_flag_table[ii][jj])
                        realvalued_flags_[ii].append(realvalued_flag_table[ii][jj])
                        uniq_mol_map_new_[ii].append(uniq_mol_map_new[ii][jj])
            print(bool_flags_)
            print(realvalued_flags_)
            print(uniq_mol_map_new_)

        Other Parameters
        ----------------
        trees (i.e., PCR tree)
            [['2', '5', '6', '12', '13', '15'],
             ['3', '4', '7', '9', '12'],
             ['2', '6', '9', '10', '11', '12', '13', '15'],
             ['2', '4', '6', '8', '11'],
             ['3', '4', '7', '9', '13', '14']]
        trees_np
            [['2' '5' '6' '12' '13' '15' '0' '0']
             ['3' '4' '7' '9' '12' '0' '0' '0']
             ['2' '6' '9' '10' '11' '12' '13' '15']
             ['2' '4' '6' '8' '11' '0' '0' '0']
             ['3' '4' '7' '9' '13' '14' '0' '0']]
        bool_flag_table
            [[ True False False False False False False False]
             [ True  True  True  True False False False False]
             [ True False False False False False False False]
             [ True  True False False False False False False]
             [ True  True  True  True False False False False]]
        realvalued_flag_table
            [['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
             ['3' '4' '7' '9' '-1' '-1' '-1' '-1']
             ['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
             ['2' '4' '-1' '-1' '-1' '-1' '-1' '-1']
             ['3' '4' '7' '9' '-1' '-1' '-1' '-1']]

        Returns
        -------

        """
        umi_map = pfreader().generic(self.args['init_seq_setting']['umi_lib_fpn'])[0].to_dict()
        # print(umi_map)

        nn = pcr_dict['data'].shape[0]
        spl_ids = rannum().uniform(
            low=0, high=nn, num=500, use_seed=False, seed=1
        )
        # print(spl_ids)
        spl_data = pcr_dict['data'][spl_ids]
        spl_id_map = self.tactic6(spl_data)
        trees = spl_data[:, 0].ravel().tolist()
        # print(trees)

        res_data = []
        mol_map_to_all_its_pcr_trees = {}
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]] = []
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]].append(tree.split('_')[1:])
        # print(mol_map_to_all_its_pcr_trees)
        for k, trees in mol_map_to_all_its_pcr_trees.items():
            print(k, trees)
            read = umi_map[int(k)]
            read_cache_table = [[] for _ in range(len(trees))]
            bool_flag_table = [[True] for _ in range(len(trees))]
            realvalued_flag_table = [[] for _ in range(len(trees))]
            trees_ori = [[] for _ in range(len(trees))]
            max_len_pcr_tree = max([len(tree) for tree in trees])
            for _, tree in enumerate(trees):
                if len(tree) < max_len_pcr_tree:
                    tree += ['0' for _ in range(max_len_pcr_tree - len(tree))]
                # print(k, tree)
            trees_np = np.array(trees)

            # ### /*** block. construct bool and real-valued flag tables ***/
            #         trees (i.e., PCR tree)
            #             [['2', '5', '6', '12', '13', '15'],
            #              ['3', '4', '7', '9', '12'],
            #              ['2', '6', '9', '10', '11', '12', '13', '15'],
            #              ['2', '4', '6', '8', '11'],
            #              ['3', '4', '7', '9', '13', '14']]
            #         trees_np
            #             [['2' '5' '6' '12' '13' '15' '0' '0']
            #              ['3' '4' '7' '9' '12' '0' '0' '0']
            #              ['2' '6' '9' '10' '11' '12' '13' '15']
            #              ['2' '4' '6' '8' '11' '0' '0' '0']
            #              ['3' '4' '7' '9' '13' '14' '0' '0']]
            #         bool_flag_table
            #             [[ True False False False False False False False]
            #              [ True  True  True  True False False False False]
            #              [ True False False False False False False False]
            #              [ True  True False False False False False False]
            #              [ True  True  True  True False False False False]]
            #         realvalued_flag_table
            #             [['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
            #              ['3' '4' '7' '9' '-1' '-1' '-1' '-1']
            #              ['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
            #              ['2' '4' '-1' '-1' '-1' '-1' '-1' '-1']
            #              ['3' '4' '7' '9' '-1' '-1' '-1' '-1']]
            # ### /*** block. construct bool and real-valued flag tables ***/
            for id_horiz_in_a_tree in range(trees_np.shape[1]):
                repeat_in_a_col = self.findListDuplicates(trees_np[:, id_horiz_in_a_tree])
                # print('repeated in a col', repeat_in_a_col)
                for id_vert_across_trees, ele_in_a_col in enumerate(trees_np[:, id_horiz_in_a_tree]):
                    if ele_in_a_col in repeat_in_a_col:
                        ids_repeat_in_a_col = [i for i, value in
                                               enumerate(trees_np[:, id_horiz_in_a_tree]) if
                                               value == ele_in_a_col]
                        # print(ids_repeat_in_a_col)
                        if bool_flag_table[id_vert_across_trees][id_horiz_in_a_tree] is False:
                            # print('repeated ele: {}'.format(ele_in_a_col))
                            bool_flag_table[id_vert_across_trees].append(False)
                            realvalued_flag_table[id_vert_across_trees].append(-1)
                        else:
                            # print('non-repeated ele: {}'.format(ele_in_a_col))
                            inspector_flags = [1 if bool_flag_table[i][id_horiz_in_a_tree] is True else 0 for i in
                                               ids_repeat_in_a_col]
                            # print(inspector_flags)
                            if sum(inspector_flags) > 1:
                                bool_flag_table[id_vert_across_trees].append(True)
                                realvalued_flag_table[id_vert_across_trees].append(ele_in_a_col)
                            else:
                                bool_flag_table[id_vert_across_trees].append(False)
                                realvalued_flag_table[id_vert_across_trees].append(-1)
                    else:
                        bool_flag_table[id_vert_across_trees].append(False)
                        realvalued_flag_table[id_vert_across_trees].append(-1)

            bool_flag_table = np.array(bool_flag_table)[:, 1:]
            realvalued_flag_table = np.array(realvalued_flag_table)
            print(trees_np)
            print(bool_flag_table)
            print(realvalued_flag_table)

            for jj in range(trees_np.shape[1]):
                read_for_repeat_tmp_per_col = {}
                for ii, val_in_a_col in enumerate(trees_np[:, jj]):
                    if val_in_a_col != '0':
                        if jj == 0:
                            # print(val_in_a_col, bool_flag_table[ii][jj])
                            if bool_flag_table[ii][jj] == True:
                                if val_in_a_col not in [*read_for_repeat_tmp_per_col.keys()]:
                                    r1 = self.mutated(
                                        read=read,
                                        pcr_error=self.args['pcr_error'],
                                    )
                                    read_for_repeat_tmp_per_col[val_in_a_col] = r1
                                    read_cache_table[ii].append(r1)
                                else:
                                    read_cache_table[ii].append(read_for_repeat_tmp_per_col[val_in_a_col])
                            else:
                                r1 = self.mutated(
                                    read=read,
                                    pcr_error=self.args['pcr_error'],
                                )
                                read_cache_table[ii].append(r1)
                        if jj > 0:
                            if bool_flag_table[ii][jj] == True:
                                if val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj])) not in [*read_for_repeat_tmp_per_col.keys()]:
                                    r1 = self.mutated(
                                        read=read_cache_table[ii][jj - 1],
                                        pcr_error=self.args['pcr_error'],
                                    )
                                    read_for_repeat_tmp_per_col[val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj]))] = r1
                                    read_cache_table[ii].append(r1)
                                else:
                                    read_cache_table[ii].append(read_for_repeat_tmp_per_col[val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj]))])
                                # print('id', [i for i, value in enumerate(trees_np[:, jj]) if value == val_in_a_col])
                            else:
                                r1 = self.mutated(
                                    read=read_cache_table[ii][jj - 1],
                                    pcr_error=self.args['pcr_error'],
                                )
                                read_cache_table[ii].append(r1)
            print(read_cache_table)

            for ii, u in enumerate(trees_np):
                for jj, v in enumerate(u):
                    if v != '0':
                        trees_ori[ii].append(trees_np[ii][jj])
            # print(trees_ori)

            for i, tree in enumerate(trees_ori):
                # print(read_cache_table[i])
                # print(read_cache_table[i][-1])
                res_data.append([
                    read_cache_table[i][-1] if read_cache_table[i] != [] else read,  # read
                    str(k) + '_' + '_'.join(tree) if read_cache_table[i] != [] else str(k),  # sam id
                    spl_id_map[str(k) + '_' + '_'.join(tree)] if read_cache_table[i] != [] else 'init',  # source
                ])
        # print(res_data)
        # print(len(res_data))
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

    def mutated(self, read, pcr_error):
        num_err_per_read = rannum().binomial(
            n=len(read), p=pcr_error, use_seed=False, seed=False
        )
        pos_list = rannum().uniform(
            low=0, high=len(read), num=num_err_per_read, use_seed=False, seed=False
        )
        base_list = rannum().uniform(
            low=0, high=3, num=num_err_per_read, use_seed=False
        )
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