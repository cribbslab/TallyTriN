__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from src.util.sequence.fastq.Read import read as rfastq
from src.util.sequence.fastq.Write import write as wfastq
from Path import to
import pandas as pd
from collections import Counter
from src.sequencing.reads.umi.Filter import filter
from src.sequencing.reads.umi.RuleOut import ruleOut as umiro
from src.sequencing.reads.seq.RuleOut import ruleOut as seqro


class template(object):

    def __init__(self, *args):
        self.args = args[0]
        self.filter = filter()
        self.umiro = umiro(read_summary=self.args)
        self.seqro = seqro(read_summary=self.args)
        self.rfastq = rfastq
        self.wfastq = wfastq

    def todf(self):
        print('reading from fastq...')
        names, seqs, _, _ = self.rfastq().fromgz(
            fastq_path=self.args['fastq']['path'],
            fastq_name=self.args['fastq']['name'],
            method='pyfastx',
        )
        df = pd.DataFrame(seqs, columns=['seq_raw'])
        df['name'] = names
        print(df)
        compo_struct = self.args['seq_struct'].split('+')
        print(compo_struct)
        umi_pos_in_struct = [i for i, cstruct in enumerate(compo_struct) if 'umi' in cstruct.split('_')]
        seq_pos_in_struct = [i for i, cstruct in enumerate(compo_struct) if 'seq' in cstruct.split('_')]
        print(umi_pos_in_struct)
        print(seq_pos_in_struct)
        umi_rule_out_len_dict = self.umiro.sequential(
            umi_pos_in_struct=umi_pos_in_struct,
            compo_struct=compo_struct,
        )
        print(umi_rule_out_len_dict)
        seq_rule_out_len_dict = self.seqro.sequential(
            seq_pos_in_struct=seq_pos_in_struct,
            compo_struct=compo_struct,
        )
        print(seq_rule_out_len_dict)
        for key, val in umi_rule_out_len_dict.items():
            start = val
            end = val + self.args[key]['len']
            df[key] = df.apply(lambda x: self.filter.singleStart(x['seq_raw'], start, end), axis=1)
        for key, val in seq_rule_out_len_dict.items():
            start = val
            end = val + self.args[key]['len']
            df[key] = df.apply(lambda x: self.filter.singleStart(x['seq_raw'], start, end), axis=1)
        print(df)
        return df

    def togz(self, df):
        if 'seq_1' not in df.columns.tolist():
            df['seq_1'] = 'B'
        df_col_names = df.columns.tolist()
        df_col_names.remove('seq_raw')
        df_col_names.remove('seq_1')
        df_trimmed = df[['seq_1'] + df_col_names]
        print(df_trimmed)
        self.wfastq().togz(
            list_2d=df_trimmed.values.tolist(),
            sv_fp=self.args['fastq']['trimmed_path'],
            fn=self.args['fastq']['trimmed_name'],
        )
        return df_trimmed


if __name__ == "__main__":
    DEFINE = {
        'umi_1': {
            'len': 12,
        },
        'umi_2': {
            'len': 10,
        },
        'umi_3': {
            'len': 12,
        },
        'seq_struct': 'umi_1+seq_1',
        # 'seq_struct': 'primer+umi_start+seq+umi_1+primer',
        # 'seq_struct': 'primer+umi_1+seq_1+primer+umi_2+seq_2+umi_3+primer',
        'primer': {
            'len': 20,
        },
        'seq_1': {
            'len': 75,
        },
        'seq_2': {
            'len': 75,
        },
        'fastq': {
            'path': to('data/'),
            'name': 'simu',
            # 'name': 'simu_u100_pcr3',
            'trimmed_path': to('data/'),
            'trimmed_name': 'simu_trimmed',
        },
        'method': 'single_start',
    }
    p = template(DEFINE)

    df = p.todf()

    print(p.togz(df))