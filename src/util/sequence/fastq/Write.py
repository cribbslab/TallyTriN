__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

import os, sys
dis = os.path.abspath("../../../")
print(dis)
sys.path.append(dis)
from Path import to
import gzip


class write(object):

    def __init__(self):
        pass

    def togz(self, list_2d, sv_fp, fn):
        fname = sv_fp + fn + '.fastq.gz'
        f = gzip.open(fname, 'wt')
        for i, read in enumerate(list_2d):
            seq = str(read[0])
            # print('No.{} saving in FASTQ format.'.format(i + 1))
            f.write('@' + '_'.join(read[1:]) + '\n')
            f.write(seq + '\n')
            f.write('+' + '\n')
            f.write('s' * len(seq) + '\n')
        f.close()
        return 0