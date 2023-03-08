import os
import pysam
import pandas as pd
import numpy as np
import itertools
import click
import collections

@click.group()
def cli():
    pass

def collapse_umi(string):
    dimers = list()
    umis = set()
    while len(string) != 0:
        try:
            dimers.append(set(string[0:3]))
            string = string[3:]
        except:
            print("umi existing indel or wrong")

    for val in itertools.product(*dimers):
        collapse = ''.join(val)
        umis.add(collapse)
    return umis

def allbasesSame(s):
    if len(set(s)) < 3:
        return True
    else:    
        return False

def collapse_cmi(string):
    dimers = list()
    umis = set()
    while len(string) != 0:
        try:
            if allbasesSame(string[0:3]):
                base = collections.Counter(string[0:3]).most_common(1)[0][0]
                dimers.append(base)
                string = string[3:]
            else:
                dimers.append(set(string[0:3]))
                string = string[3:]
        except:
            print("umi existing indel or wrong")

    for val in itertools.product(*dimers):
        collapse = ''.join(val)
        umis.add(collapse)
    return umis

def umi_greedy(input_umis): 
    merged_umis = {}
    merged_umis_idx = {}
    n_merged = len(input_umis)
    left_umis = input_umis            
    n_steps=0
    while n_merged > 1:
        umi_ids = dict()
        for ii, uu_list in enumerate(left_umis):
            for uu in uu_list:
                if uu in umi_ids:
                    umi_ids[uu].append(ii)
                else:
                    umi_ids[uu] = [ii]
        umi_counts = {kk:len(vv) for kk, vv in umi_ids.items()}
        umi_df = pd.DataFrame(
            {'umi_id':umi_counts.keys(), 'umi_count':umi_counts.values()}
        ).sort_values('umi_count', ascending=False).reset_index(drop=True)

        n_merged = umi_df.umi_count[0]
        if n_merged > 1:
            merged_umis_idx[n_steps] = umi_ids[umi_df.umi_id[0]]
            merged_umis[n_steps] = umi_df.umi_id[0]
        # fix umis in path
            ll_umis = []
            for rr_idx, ll_uu in enumerate(left_umis):
                if rr_idx not in merged_umis_idx[n_steps]:
                    ll_umis.append(ll_uu)
            left_umis = ll_umis
            if len(left_umis) == 0:
                break
            n_steps += 1
        else:
            break
    shortest_path = len(input_umis) - sum([len(ii) - 1 for ii in merged_umis_idx.values()])
    return(shortest_path)

@cli.command()
@click.option("-i", "--inbam", help="the input bam file")
@click.option("-t", "--tag", help="The tag add to genes or transcripts")
@click.option("-o", "--output", help="The output file for genes or transcripts-cell count")
def count(inbam, tag, output, sep = "_"):
    """Here read the bam file and extract information of cell barcodes, umis, and transcripts with tags"""
    tab = dict()
    count = 0
    n_tag = 0
    outf = open(output, "w")
    outf.write('%s\t%s\t%s\n' % ('gene', 'cell', 'count'))
    with pysam.AlignmentFile(inbam) as bf:
        for i, r in enumerate(bf):
            if r.has_tag(tag) is True:
                key = (r.get_tag(tag), r.qname.split(sep)[1])
                tab.setdefault(key,[]).append(r.qname.split(sep)[2])
                n_tag += 1
            else:
                pass
        
        print("The total number of input reads is ", i+1)
        print("The total number of input reads with XT tag is ", n_tag)
        
        n_reads = 0
        n = 0
        ii = 0
        m = 0     # how many expressed genes have only one read
        mst = 0   #how many expressed genes are duplicated into 1 by majority rule
        maj_cr = 0     # how many reads are corrected by majority rule
        ilp_cr = 0  # how many reads are corrected by majority rule and greedy
        for kk, vv in tab.items():
            if not ii % 3000:
                print(ii)
            if len(vv) == 1:
                count = 1
                n += count
                m += 1
                outf.write('%s\t%s\t%s\n' % (kk[0], kk[1], count))
            else:
                #corrected_cmis = [collapse_cmi(uu) for uu in vv]
                corrected_cmis = set(tuple(collapse_cmi(uu)) for uu in vv)
                    #use_ilp = any([len(xx) > 1 for xx in corrected_cmis])                
                    #if not use_ilp:
                if len(corrected_cmis) == 1:
                    count = 1
                    n += count
                    mst += 1
                    maj_cr += len(vv)
                else:
                    count = umi_greedy([collapse_umi(uu) for uu in vv])
                    if count == 1:
                        n_reads += len(vv)
                    n += count
                    ilp_cr += len(vv)
                outf.write('%s\t%s\t%s\n' % (kk[0], kk[1], count))
            ii += 1
    print("The total count of umi_count is ", n)
    print("The number of gene expression itself is only 1 is: ", m)
    print("Of all gene expression, the percentage of gene expression itself that is only 1 is %d / %d = " % (m, ii, ), m/(ii))
    print("The number of gene expression quantified as 1 by majority rule is: ", mst)
    print("How many reads are corrected only by majority: ", maj_cr)
    print("How many reads are corrected by greedy after majority: ", ilp_cr)
    print("The total number of reads corresponding to gene expression quantified as 1 by greedy is: ", n_reads)
    return

if __name__ == "__main__":
    cli()
