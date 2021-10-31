import regex


class patternMatching(object):

    def __init__(self, ):
        pass

    def regexCompile(self, cond='{e<=1}', query='TTTTTTTTRTTTTTTTTTTT', std_pattern='TTTTTTTTTTTTTTTTTTTT'):
        r = regex.compile(''.join(['(%s)', cond]) % std_pattern)
        res = r.match(query)
        if res:
            return True
        else:
            return False

    def regexfindall(self, cond='{e<=1}', query='TTTTTTTTRTTTTTTTTTTT', std_pattern='TTTTTTTTTTTTTTTTTTTT'):
        res = regex.findall(''.join(['(', std_pattern, ')', cond]), query)
        # print(res)
        return len(res)
        # if len(res):
        #     return True
        # else:
        #     return False


if __name__ == "__main__":
    p = patternMatching()

    # print(p.regexCompile(
    #     query='AAGCAGTGGTATCAACACAGAAATTACTTTCCCAAACCCCCCAAAGGAAAGGAAACCCGGGTTTTTTTTTTTTTTTTTTTTTTTTTTGGGAAAATTCACTCTGCGTTGATACCACTGCTAATTAAGGGGGATTCACTCTGCGTTGATACCACTGCTTGGTTTCAAGGGGGATTCACTCTGCGTTGATACCACTGCTTGGTTTCCAGGGGGATTCACTCTGCGTTGATACCACTGCTTAGCAATGCATGTG',
    #     std_pattern="TTTTTTTTTTTTTTTTTTTT"
    # ))

    print(p.regexfindall(
        query='ATCGATCGGCATGCAGTGCAGAACTGACGAT',
        std_pattern="ATCT",
        cond='{s<=0}'
    ))