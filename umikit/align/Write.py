__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

# import pysam
from umikit.util.Console import console


class write(object):

    def __init__(self, df, is_sv=False, verbose=False):
        self.df = df
        self._is_sv = is_sv
        self.console = console()
        self.console.verbose = verbose

    @property
    def is_sv(self, ):
        return self._is_sv

    @is_sv.setter
    def is_sv(self, value):
        self._is_sv = value

    def tobam(self, tobam_fpn, tmpl_bam_fpn, whitelist=[]):
        if self._is_sv:
            tmpl_bam = pysam.AlignmentFile(tmpl_bam_fpn, "rb")
            write_to_bam = pysam.AlignmentFile(tobam_fpn, "wb", template=tmpl_bam)
            fs = self.df.loc[self.df['id'].isin(whitelist)]['read']
            for i in fs:
                # print(i)
                write_to_bam.write(i)
            write_to_bam.close()
            return write_to_bam
        else:
            return 'BAM is not saved because of no assignment from the user.'