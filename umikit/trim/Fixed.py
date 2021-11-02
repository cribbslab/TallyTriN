__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import argparse
from umikit.trim.Template import template as umitrim
from umikit.util.Console import console


class fixed():

    def __init__(self, mode='external', params=None, verbose=True):
        if mode == 'internal':
            self.params = params
            self.verbose = verbose
            self.console = console()
            self.console.verbose = self.verbose
            self.console.print('run Mclumi internally.')
            self.console.print('Your params for trimming UMIs are: \n{}'.format(self.params))
        else:
            self.parser = argparse.ArgumentParser(description='Sequence Identity calculations')
            self.parser.add_argument(
                "--read_structure", "-rs",
                metavar='read_structure',
                dest='rs',
                required=True,
                type=str,
                help='the read structure with elements in conjunction with +',
            )
            self.parser.add_argument(
                "--lens", "-l",
                metavar='lens',
                dest='l',
                required=True,
                type=str,
                help='lengths of all sub-structures separated by +',
            )
            self.parser.add_argument(
                "--input", "-i",
                metavar='input',
                dest='i',
                required=True,
                type=str,
                help='input a fastq file in gz format for trimming UMIs',
            )
            self.parser.add_argument(
                "--output", "-o",
                metavar='output',
                dest='o',
                required=True,
                type=str,
                help='output a UMI-trimmed fastq file in gz format.',
            )
            self.parser.add_argument(
                "--verbose", "-vb",
                metavar='verbose',
                dest='vb',
                default=True,
                type=bool,
                help='to enable if output logs are on console',
            )
            args = self.parser.parse_args()
            self.verbose = args.vb
            self.read_structure = args.rs
            self.structure_lengths = []
            for i in args.l.split('+'):
                self.structure_lengths.append(int(i))
            # print(args.l)
            # print(args.i)
            self.fastq_fpn = args.i
            self.fastq_trimmed_fpn = args.o
            self.params = {s: {'len': l} for s, l in zip(self.read_structure.split('+'), self.structure_lengths)}
            self.params['read_struct'] = self.read_structure
            self.params['fastq'] = {
                'fpn': self.fastq_fpn,
                'trimmed_fpn': self.fastq_trimmed_fpn,
            }
            self.console = console()
            self.console.verbose = self.verbose

    def call(self, ):
        umitrim_parser = umitrim(self.params, self.verbose)
        df = umitrim_parser.todf()
        umitrim_parser.togz(df)


if __name__ == "__main__":
    from Path import to

    params = {
        'umi_1': {
            'len': 12,
        },
        'umi_2': {
            'len': 10,
        },
        'umi_3': {
            'len': 12,
        },
        # 'read_struct': 'umi_1+seq_1',
        'read_struct': 'primer_1+umi_1+seq_1+seq_2+umi_2+umi_3+primer_2',
        # 'read_struct': 'primer+umi_start+seq+umi_1+primer',
        # 'read_struct': 'primer+umi_1+seq_1+primer+umi_2+seq_2+umi_3+primer',
        'primer_1': {
            'len': 20,
        },
        'primer_2': {
            'len': 20,
        },
        'seq_1': {
            'len': 6,
        },
        'seq_2': {
            'len': 8,
        },
        'fastq': {
            'fpn': to('example/data/pcr_1.fastq.gz'),
            'trimmed_fpn': to('example/data/pcr_1_trim.fastq.gz'),
        },
    }
    p = fixed(
        mode='internal',
        # mode='external',

        params=params,
    )
    p.call()