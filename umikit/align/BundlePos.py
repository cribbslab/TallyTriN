# import umi_tools.sam_methods as sam_methods
import mclumi.external.UMItoolsSamMethods as sam_methods
import pysam
from umikit.util.Hamming import hamming
import pandas as pd
import numpy as np
from Path import to


def edave(x, d):
    repr_nodes = d[x[1]]
    node_len = len(repr_nodes)
    if node_len != 1:
        ed_list = []
        for i in range(node_len):
            for j in range(i + 1, node_len):
                ed_list.append(hamming().general(
                    s1=repr_nodes[i],
                    s2=repr_nodes[j]
                ))
        return np.ceil(sum(ed_list) / (len(ed_list)))
    else:
        return -1


def convert(options, in_fpn, out_fpn):
    bundle_iterator = sam_methods.get_bundles(
        options,
        metacontig_contig=None,
    )
    infile = pysam.Samfile(in_fpn, 'rb')
    inreads = infile.fetch()
    num_bundled_pos = 0
    umis = []
    nInput = 0
    reads = []
    uniq_umi_cnt = 0
    write_to_bam = pysam.AlignmentFile(out_fpn, "wb", template=infile)
    for i, (bundle, key, status) in enumerate(bundle_iterator(inreads)):
        # print([bundle[umi]["read"] for umi in bundle])
        for j, umi in enumerate(bundle):
            read = bundle[umi]["read"]
            read.set_tag('PO', i)
            for _ in range(bundle[umi]["count"]):
                reads.append(bundle[umi]["read"])
            uniq_umi_cnt += 1
        nInput += sum([bundle[umi]["count"] for umi in bundle])

        umis.append([i.decode('utf-8') for i in bundle.keys()])
        num_bundled_pos += 1

    for y in reads:
        write_to_bam.write(y)
    write_to_bam.close()

    print('# of unique reads {}'.format(len(reads)))
    print('# of bundled pos {}'.format(num_bundled_pos))
    print('# of repeated UMIs in total {}'.format(nInput))
    df = pd.DataFrame(index=np.arange(num_bundled_pos))
    df[1] = np.arange(num_bundled_pos)
    df[2] = df.apply(lambda x: edave(x, umis), axis=1)
    return df[2].value_counts()


in_fpn = to('example/data/example.bam')
out_fpn = to('example/data/example_bundle.bam')
options = {'stats': 'deduplicated',
           'get_umi_method': 'read_id',
           'umi_sep': '_',
           'umi_tag': 'RX',
           'umi_tag_split': None,
           'umi_tag_delim': None,
           'cell_tag': None,
           'cell_tag_split': '-',
           'cell_tag_delim': None,
           'filter_umi': None,
           'umi_whitelist': None,
           'umi_whitelist_paired': None,
           'method': 'directional',
           'threshold': 1,
           'spliced': False,
           'soft_clip_threshold': 4,
           'read_length': False,
           'per_gene': False,
           'gene_tag': None,
           'assigned_tag': None,
           'skip_regex': '^(__|Unassigned)',
           'per_contig': False,
           'gene_transcript_map': None,
           'per_cell': False,
           'whole_contig': False,
           'detection_method': None,
           'mapping_quality': 0,
           'output_unmapped': False,
           'unmapped_reads': 'discard',
           'chimeric_pairs': 'use',
           'unpaired_reads': 'use',
           'ignore_umi': False,
           'ignore_tlen': False,
           'chrom': None,
           'subset': None,
           'in_sam': False,
           'paired': False,
           'out_sam': False,
           'no_sort_output': False,
           'stdin': "<_io.TextIOWrapper name='example.bam' mode='r' encoding='UTF-8'>",
           'stdlog': "<_io.TextIOWrapper name='<stdout>' mode='w' encoding='UTF-8'>", 'log2stderr': False,
           'compresslevel': 6,
           'timeit_file': None,
           'timeit_name': 'all',
           'timeit_header': None,
           'loglevel': 1,
           'short_help': None,
           'random_seed': None
           }
convert(options=options, in_fpn=in_fpn, out_fpn=out_fpn)