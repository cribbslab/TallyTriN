__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import pyfastx


class read():

    def __init__(self):
        pass

    def fromgz(self, fastq_fpn):
        """

        Note
        ----
        read and parse a Fastq file.

        Parameters
        ----------
        fastq_fpn
            the path of a fastq file
        Returns
        -------
            tuple consisting of names, seqs, placeholders, qualities
        """
        names = []
        seqs = []
        placeholders = []
        qualities = []
        fq = pyfastx.Fastx(fastq_fpn)
        for name, seq, qual, comment in fq:
            seqs.append(''.join(seq))
            names.append(''.join(name))
        return names, seqs, placeholders, qualities
