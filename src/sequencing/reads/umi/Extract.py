__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

from src.util.sequence.fastq.Read import read as rfastq
from Path import to
import pandas as pd
from collections import Counter


class extract(object):

    def __init__(self, *args):
        self.args = args[0]
        self.rfastq = rfastq

    def filter(self, x, rule):
        a_len = 0
        b_len = 0
        for j in rule[0]:
            b_len += self.args[j]['len']
        for j in rule[1]:
            a_len += self.args[j]['len']
        return x[b_len: b_len+self.args['umi']['len']]
        # return x[a_len: a_len+self.args['umi']['len']]

    def cus(self):
        print('reading...')
        seqs = self.rfastq().fromgz(
            fastq_path=self.args['fastq']['path'],
            fastq_name=self.args['fastq']['name'],
        )
        df_seq = pd.DataFrame(seqs)
        dd = self.args['seq_struct'].split('*')
        umi_pos = [i for i, d in enumerate(dd) if d == 'umi']
        rule = []
        if len(umi_pos) != 1:
            rule.append(dd[:umi_pos[0]])
            rule.append(dd[umi_pos[1]+1:])
        else:
            rule.append(dd[:umi_pos[0]])
        umis = df_seq.apply(lambda x: self.filter(x[0], rule), axis=1)
        return umis


if __name__ == "__main__":
    DEFINE = {
        'umi': {
            'len': 36,
        },
        'seq_struct': 'primer*umi*seq*umi*primer',
        'primer': {
            'len': 20,
        },
        'fastq': {
            'path': to('data/'),
            'name': 'simu',
        },
    }
    p = extract(DEFINE)

    umis = p.cus()

    print()