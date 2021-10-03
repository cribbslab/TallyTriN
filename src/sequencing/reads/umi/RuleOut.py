__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from Path import to


class ruleOut(object):

    def __init__(self, read_summary):
        self.read_summary = read_summary

    def sequential(self, umi_pos_in_struct, compo_struct):
        """
            if i == 0:
                rule_out_struct_dict['umi_1'] = compo_struct[:umi_pos_in_struct[i]]
            # elif i == len(umi_pos_in_struct) - 1:
            #     rule_out_struct_dict['umi_end'] = compo_struct[umi_pos_in_struct[i] + 1:]
            else:
                rule_out_struct_dict['umi_' + str(i)] = compo_struct[umi_pos_in_struct[i - 1] + 1: umi_pos_in_struct[i]]

        :param umi_pos_in_struct:
        :param compo_struct:
        :return:
        """
        rule_out_struct_dict = {}
        for i in range(len(umi_pos_in_struct)):
            if i == 0:
                rule_out_struct_dict['umi_' + str(i + 1)] = compo_struct[:umi_pos_in_struct[i]]
            else:
                rule_out_struct_dict['umi_' + str(i + 1)] = compo_struct[umi_pos_in_struct[i - 1] + 1: umi_pos_in_struct[i]]
        print('sss', rule_out_struct_dict)
        rule_out_rel_len_dict = {}
        for key, val in rule_out_struct_dict.items():
            rule_out_rel_len_dict[key] = 0
            if val == []:
                rule_out_rel_len_dict[key] = 0
            else:
                for j in val:
                    rule_out_rel_len_dict[key] += self.read_summary[j]['len']
        print(rule_out_rel_len_dict)
        rule_out_accumu_len_dict = {}
        accumu = []
        for key, val in rule_out_rel_len_dict.items():
            print(key)
            print(self.read_summary)
            accumu.append(val)
            rule_out_accumu_len_dict[key] = sum(accumu)
            accumu.append(self.read_summary[key]['len'])
        print(rule_out_accumu_len_dict)
        return rule_out_accumu_len_dict


if __name__ == "__main__":
    DEFINE = {
        'umi': {
            'len': 12,
        },
        # 'seq_struct': 'umi*seq',
        'seq_struct': 'primer*umi*seq*umi*primer',
        'primer': {
            'len': 20,
        },
        'seq': {
            'len': 75,
        },
        'fastq': {
            'path': to('data/'),
            'name': 'simu',
        },
    }
    p = ruleOut(DEFINE)

    umis = p.cus()

    print()