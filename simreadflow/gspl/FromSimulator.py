__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import rpy2.robjects as rob
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


class fromSimulator(object):

    def __init__(self, simulator):
        self.simulator = simulator

    def SPsimSeqFixSM(self, ):
        res = rob.r(
            '''
            suppressPackageStartupMessages(library(SPsimSeq))
            cat("SPsimSeq package version", as.character(packageVersion("SPsimSeq")), "\n")

            # load the Zhang bulk RNA-seq data
            data("zhang.data.sub")
            # filter genes with sufficient expression (important step) 
            zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]

            set.seed(6452)
            zhang.counts2 <- zhang.counts[sample(nrow(zhang.counts), 20), ]
            sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2,
                                      group = zhang.data.sub$MYCN.status, n.genes = 20, 
                                      batch.config = 1,
                                      group.config = c(0.5, 0.5), tot.samples = 2,
                                      pDE = 0.5, lfc.thrld = 0.5, 
                                      result.format = "list")
            sim.data.bulk1 <- sim.data.bulk[[1]]                              
            return (data.frame(sim.data.bulk1$counts))
            '''
        )
        with localconverter(rob.default_converter + pandas2ri.converter):
            df = rob.conversion.rpy2py(res)
        return df.T

    def tool(self, ):
        return {
            'SPsimSeqFixSM': self.SPsimSeqFixSM()
        }

    def run(self, ):
        return self.tool()[self.simulator]