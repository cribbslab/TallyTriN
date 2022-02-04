__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from umikit.util.Console import console


class seqRuleOut(object):

    def __init__(self, read_summary, verbose=True):
        self.read_summary = read_summary
        self.verbose = verbose
        self.console = console()
        self.console.verbose = self.verbose

    def sequential(self, compo_struct, seq_pos_in_struct):
        """

        Notes
        -----------
        Starting positions of all genomic sequences.

        Example
        -------
        rule_out_struct_dict returns all structures before each key of rule_out_struct_dict: {
            'seq_1': [struct_1, struct_2, ..., struct_n],
            'seq_2': [struct_1, struct_2, ..., struct_n],
            ...,
            'seq_m': [struct_1, struct_2, ..., struct_n],
            }
            e.g., {'seq_1': ['primer_1', 'umi_1'], 'seq_2': []}
            if the structure is 'primer_1+umi_1+seq_1+seq_2+umi_2+primer_2'
            * each key does not count all structures of its preceding keys.

        rule_out_rel_len_dict returns accumulated lengths of all structures in the list w.r.t each
        key of rule_out_struct_dict.

        rule_out_accumu_len_dict returns the starting positions of all UMIs.

        Parameters
        ----------
        compo_struct
            1d list of strings of seq_struct split by +.
        seq_pos_in_struct
            1d list of indices of compo_struct.

        Returns
        -------
        1d dict: {umi_1: int, umi_2: int, ..., umi_n: int}

        """
        self.console.print('======>finding the starting positions of all genomic sequence...')
        rule_out_struct_dict = {}
        for i in range(len(seq_pos_in_struct)):
            if i == 0:
                rule_out_struct_dict['seq_' + str(i + 1)] = compo_struct[:seq_pos_in_struct[i]]
            else:
                rule_out_struct_dict['seq_' + str(i + 1)] = compo_struct[seq_pos_in_struct[i - 1] + 1: seq_pos_in_struct[i]]
        rule_out_rel_len_dict = {}
        for key, val in rule_out_struct_dict.items():
            rule_out_rel_len_dict[key] = 0
            if val == []:
                rule_out_rel_len_dict[key] = 0
            else:
                for j in val:
                    rule_out_rel_len_dict[key] += self.read_summary[j]['len']
        rule_out_accumu_len_dict = {}
        accumu = []
        for key, val in rule_out_rel_len_dict.items():
            accumu.append(val)
            rule_out_accumu_len_dict[key] = sum(accumu)
            accumu.append(self.read_summary[key]['len'])
        for k, v in rule_out_accumu_len_dict.items():
            self.console.print('=========>{} starting position: {}'.format(k, v))
        return rule_out_accumu_len_dict