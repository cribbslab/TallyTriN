__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"


class umi(object):

    def __init__(self, ):
        self.umi_unit_len_fixed = 12
        self.umi_unit_pattern = 1
        self.pcr_errs, self.seq_errs = self.errors()
        print(self.pcr_errs)
        print(self.seq_errs)

    def errors(self, ):
        pcr_errs = []
        seq_errs = []
        e = 1e-5
        while e < 3e-1:
            pcr_errs.append(e)
            seq_errs.append(e)
            if 5 * e < 3e-1:
                pcr_errs.append(2.5 * e)
                pcr_errs.append(5 * e)
                pcr_errs.append(7.5 * e)
                seq_errs.append(2.5 * e)
                seq_errs.append(5 * e)

                seq_errs.append(7.5 * e)
            e = 10 * e
        pcr_errs.append(0.2)
        seq_errs.append(0.2)
        pcr_errs.append(0.3)
        seq_errs.append(0.3)
        # print(pcr_errs)
        # print(seq_errs)
        return pcr_errs, seq_errs

    def pcrNums(self, ):
        pass

    def pcrErrs(self, ):
        pass

    def seqErrs(self, ):
        pass

    def umiLens(self, ):
        pass

    def amplRates(self, ):
        pass


if __name__ == "__main__":
    p = umi()

    # print(p.pcrNums())

    # print(p.pcrErrs())

    # print(p.seqErrs())

    # print(p.umiLens())

    # print(p.amplRates())