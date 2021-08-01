__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__author__ = "Adam Cribbs lab"

from src.sequencing.reads.umi.Extract import extract as umiextr
from Path import to
from src.util.file.read.Reader import reader as gfreader
from collections import Counter


class reoccurence(object):

    def __init__(self, *args):
        self.args = args[0]
        self.umiextr = umiextr
        self.gfreader = gfreader()

    def pinrun(self, ):
        umis = self.umiextr(self.args).cus()
        # print(umis.unique().tolist())
        umi_lib = self.gfreader.generic(df_fpn=self.args['umi']['lib_path'])[0].tolist()
        umi_diff = list(set(umis).difference(umi_lib))
        umi_diff_reduced = self.find(umi_diff)
        # print(umi_diff_reduced)
        umi_diff_check = list(set(umi_diff_reduced).difference(umi_lib))
        # print(umi_diff_check)
        return umi_diff_check

    def find(self, arr):
        p = [i for i in range(36) if i%3 == 0]
        print(p)
        arr_ = []
        for i, m in enumerate(arr):
            arr_.append(list(m))
        for i, m in enumerate(arr):
            for t in p:
                y = Counter(m[t: t+3])
                if len(y) != 1:
                    key = list(y.keys())
                    if 2 in list(y.values()):
                        id_two = list(y.values()).index(2)
                        id_one = list(y.values()).index(1)
                        na_two = key[id_two]
                        na_one = key[id_one]
                        x = m[t: t+3].index(na_one)
                        pos = t + x
                        arr_[i][pos] = na_two
                    else:
                        continue
        arr__ = []
        for i, m in enumerate(arr_):
            arr__.append(''.join(m))
        return arr__


if __name__ == "__main__":
    DEFINE = {
        'umi': {
            'len': 36,
            'lib_path': to('data/umi.txt')
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
    p = reoccurence(DEFINE)
    print(p.pinrun())