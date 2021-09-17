__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import gzip
import os, sys
import pandas as pd
dis = os.path.abspath("../../../")
sys.path.append(dis)
from functools import wraps
from src.util.file.create.Folder import folder as crtfolder


class write(object):

    def __init__(self):
        self.crtfolder = crtfolder()

    def __call__(self, deal):
        generic = self.generic
        excel = self.excel
        @wraps(deal)
        def write(ph, *args, **kwargs):
            res = deal(ph, **kwargs)
            keys = [*kwargs.keys()]
            if kwargs['type'] == 'generic':
                generic(
                    df=kwargs['df'],
                    sv_fpn=kwargs['sv_fpn'],
                    df_sep='\t' if 'df_sep' not in keys else kwargs['df_sep'],
                    id_from=0 if 'id_from' not in keys else kwargs['id_from'],
                    header=None if 'header' not in keys else kwargs['header'],
                    index=False if 'index' not in keys else kwargs['index'],
                )
            elif kwargs['type'] == 'excel':
                excel(
                    df=kwargs['df'],
                    sv_fpn=kwargs['sv_fpn'],
                    sheet_name='Sheet1' if 'sheet_name' not in keys else kwargs['sheet_name'],
                    id_from=0 if 'id_from' not in keys else kwargs['id_from'],
                    header=None if 'header' not in keys else kwargs['header'],
                    index=False if 'index' not in keys else kwargs['index'],
                )
            return res
        return write

    def togz(self, list_2d, sv_fp, fn, symbol='_'):
        self.crtfolder.osmkdir(sv_fp)
        fname = sv_fp + fn + '.fastq.gz'
        f = gzip.open(fname, 'wt')
        for i, read in enumerate(list_2d):
            seq = str(read[0])
            # print('No.{} saving in FASTQ format.'.format(i + 1))
            f.write('@' + symbol.join(read[1:]) + '\n')
            f.write(seq + '\n')
            f.write('+' + '\n')
            f.write('s' * len(seq) + '\n')
        f.close()
        return 0

    def generic(self, df, sv_fpn, df_sep='\t', header=None, index=False, id_from=0):
        df_ = pd.DataFrame(df)
        df_.index = df_.index + id_from
        return df_.to_csv(
            sv_fpn,
            sep=df_sep,
            header=header,
            index=index
        )

    def excel(self, df, sv_fpn=None, sheet_name='Sheet1', header=None, index=False, id_from=0):
        df_ = pd.DataFrame(df)
        df_.index = df_.index + id_from
        return df_.to_excel(
            sv_fpn,
            sheet_name=sheet_name,
            header=header,
            index=index
        )