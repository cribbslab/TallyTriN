__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import os


class folder(object):

    def __init__(self, ):
        pass

    def osmkdir(self, DIRECTORY):
        """

        Parameters
        ----------
        DIRECTORY

        Returns
        -------

        """
        if not os.path.exists(DIRECTORY):
            os.makedirs(DIRECTORY)
        return 0
