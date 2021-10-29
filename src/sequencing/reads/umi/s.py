import pysam
from Path import to


class sdasd(object):

    def __init__(self, ):
        pass

    def asds(self, ):
        samfile = pysam.AlignmentFile(to('data/RM82CLK1_S1_XT_gene.bam'), "rb")
        for i, read in enumerate(samfile):
            # print(i)
            if read.tags[10][0] == 'XS':
                # print(read)
                print(read.get_tag('XT'))
                # break
        # print(samfile.header)
        return


if __name__ == "__main__":
    p = sdasd()

    print(p.asds())