__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from umikit.util.Console import console


class umiRuleOut(object):

    def __init__(self, read_summary, verbose=True):
        self.read_summary = read_summary
        self.verbose = verbose
        self.console = console()
        self.console.verbose = self.verbose

    def sequential(self, compo_struct, umi_pos_in_struct):
        """
        
        Note
        ----
        Starting positions of all UMIs.

        Example
        -------
        rule_out_struct_dict returns all structures before each key of rule_out_struct_dict: {
            'umi_1': [struct_1, struct_2, ..., struct_n],
            'umi_2': [struct_1, struct_2, ..., struct_n],
            ...,
            'umi_m': [struct_1, struct_2, ..., struct_n],
            }
            e.g., {'umi_1': ['primer_1'], 'umi_2': ['seq_1']}

        rule_out_rel_len_dict returns accumulated lengths of all structures in the list w.r.t each key of
        rule_out_struct_dict.

        rule_out_accumu_len_dict returns the starting positions of all UMIs.
        
        Parameters
        ----------
        compo_struct
            1d list of strings of seq_struct split by +.
        umi_pos_in_struct
            1d list of indices of compo_struct.

        Returns
        -------
        1d dict: {umi_1: int, umi_2: int, ..., umi_n: int}

        """
        self.console.print('======>finding the starting positions of all UMIs...')
        rule_out_struct_dict = {}
        for i in range(len(umi_pos_in_struct)):
            if i == 0:
                rule_out_struct_dict['umi_' + str(i + 1)] = compo_struct[:umi_pos_in_struct[i]]
            else:
                rule_out_struct_dict['umi_' + str(i + 1)] = compo_struct[umi_pos_in_struct[i - 1] + 1: umi_pos_in_struct[i]]
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