__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "GPL v3.0"
__lab__ = "Adam Cribbs lab"

from src.sequencing.reads.simulate.inf.Pseudo import pseudo as seqpseudo
from src.sequencing.reads.umi.Library import library as liblogginger


class design(seqpseudo):

    def __init__(self, *args, **kwargs):
        super(design, self).__init__(*args, **kwargs)
        self.args = args
        self.kwargs = kwargs

    @liblogginger()
    def general(self, lib_fpn='./spacer.txt', is_sv=True):
        return ''.join([
            self.kwargs['dna_map'][i] for i in
            self.kwargs['pseudorandom_num']
        ])