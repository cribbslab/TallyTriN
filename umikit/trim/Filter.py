__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from Path import to


class filter(object):

    def __init__(self, ):
        pass

    def method(self, ):
        return {
            'single_start': self.singleStart,
            'paired': self.paired,
        }

    def singleStart(self, x, start, end):
        """

        Parameters
        ----------
        x
        start
        end

        Returns
        -------

        """
        return x[start: end]
        # len_dict = {}
        # for key, val in start_len.items():
        #     len_dict[key] = 0
        #     for j in val:
        #         len_dict[key] += self.read_summary[j]['len']
        # print(len_dict)
        # # return x[b_len: b_len+self.read_summary['umi']['len']]
        # # return x[a_len: a_len+self.args['umi']['len']]

    def paired(self, x, rule):
        """

        Parameters
        ----------
        x
        rule

        Returns
        -------

        """
        a_len = 0
        b_len = 0
        for j in rule[0]:
            b_len += self.args[j]['len']
        for j in rule[1]:
            a_len += self.args[j]['len']
        return x[b_len: b_len+self.args['umi']['len']]
        # return x[a_len: a_len+self.args['umi']['len']]


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
            'len': 20,
        },
        'fastq': {
            'path': to('data/'),
            'name': 'simu',
        },
    }
    p = filter(DEFINE)

    umis = p.cus()

    print()