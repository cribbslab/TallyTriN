"""===========================
Pipeline count
===========================

Overview
========

The aim of this pipeline is to take a nanopore input fastq and then process
the file so a counts matrix is generated for downstream differential expression.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_cat_fastq.py config

Input files
-----------

fastq.gz file of nanopore reads that have been sequenced with trimers at
the 5' and 3' end. Data should be added to the data.dir folder.

Pipeline output
===============

A counts matrix with sample as columns and rows as either transcripts or genes. 


Code
====

"""
import sys
import os
import pysam
import glob
import pandas as pd
from ruffus import *
import cgatcore.iotools as iotools
import cgatcore.pipeline as P
import cgatcore.experiment as E
from cgatcore.pipeline import cluster_runnable

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


SEQUENCESUFFIXES = ("*.fastq.gz")

FASTQTARGET = tuple([os.path.join("data.dir/", suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])


def merge_feature_data(infiles):
    '''will merge all of the input files'''

    final_df = pd.DataFrame()
    for infile in infiles:
        name = infile.replace(".counts.tsv.gz", "")
        tmp_df = pd.read_table(infile, sep="\t", header=0, names=["transcript_name", name], index_col=0)
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True)

    return final_df

def merge_featurecounts_data(infiles):
    '''will merge all of the input files from featurecounts count output'''

    final_df = pd.DataFrame()
    for infile in infiles:
    
        tmp_df = pd.read_table(infile, sep="\t", header=0, index_col=0, skiprows = 1)
        tmp_df = tmp_df.iloc[:,-1:]
        tmp_df.columns = ["count"]
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True, suffixes=("","_drop"))

    names = [x.replace("_gene_assigned.txt", "") for x in infiles]
    final_df.columns = names
    return final_df

# fastqsplitter input.fastq.gz -n 3 --prefix split

@follows(mkdir("processed_fastq.dir"))
@transform(FASTQTARGET,
         regex("data.dir/(\S+).fastq.gz"),
         r"processed_fastq.dir/\1_polyA.fastq.gz")
def polya_correct(infile, outfile):
    '''filter less than 300 bp reads and then make sure polyA is in correct orientation'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/complement_polyA.py --infile=%(infile)s --outname=%(outfile)s'''

    P.run(statement)

@transform(polya_correct,
         regex("processed_fastq.dir/(\S+)_polyA.fastq.gz"),
         r"processed_fastq.dir/\1_tso_UMI.fastq.gz")
def polya_umi(infile, outfile):
    '''Identify the polya umi for each read'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    if not PARAMS['correct']:
        statement = """python %(PYTHON_ROOT)s/polya_umi_nocorrect.py --infile=%(infile)s --outname=%(outfile)s"""

    else:
        statement = '''python %(PYTHON_ROOT)s/polya_umi.py --infile=%(infile)s --outname=%(outfile)s'''

    P.run(statement)


@transform(polya_umi,
         regex("processed_fastq.dir/(\S+)_tso_UMI.fastq.gz"),
         r"processed_fastq.dir/\1_polya_tso_UMI.fastq.gz")
def tso_umi(infile, outfile):
    '''Identify the tso umi for each read'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    if PARAMS['tso_present']:

        if not PARAMS['correct']:
            statement = """python %(PYTHON_ROOT)s/tso_umi_nocorrect.py --infile=%(infile)s --outname=%(outfile)s"""

        else:
            statement = '''python %(PYTHON_ROOT)s/tso_umi.py --infile=%(infile)s --outname=%(outfile)s'''

    else:
        statement = """cp %(infile)s %(outfile)s"""

    P.run(statement)


@follows(mkdir("mapped_files.dir"))
@transform(tso_umi,
           regex("processed_fastq.dir/(\S+)_polya_tso_UMI.fastq.gz"),
           r"mapped_files.dir/\1_tso_polya_UMI.sam")
def mapping_trans(infile, outfile):
    '''map using minimap2 for the transcripts'''

    statement = '''minimap2 -ax map-ont -p 0.9 --end-bonus 10 -N 3 %(cdna_fasta)s %(infile)s  > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement)


@transform(mapping_trans,
           regex("mapped_files.dir/(\S+)_tso_polya_UMI.sam"),
           r"mapped_files.dir/\1_final_sorted.bam")
def samtools(infile, outfile):
    '''run samtools on the output and index'''

    name = outfile.replace("_final_sorted.bam", "")

    statement = '''samtools view -bS %(infile)s > %(name)s_final.bam && 
                   samtools sort %(name)s_final.bam -o %(name)s_final_sorted.bam && 
                   samtools index %(name)s_final_sorted.bam '''

    P.run(statement)


@transform(samtools,
           regex("mapped_files.dir/(\S+)_final_sorted.bam"),
           r"mapped_files.dir/\1_XT.bam")
def xt_tag(infile, outfile):
    '''add XT tag to the samfile'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/add_XT.py --infile=%(infile)s --outname=%(outfile)s && samtools index %(outfile)s'''

    P.run(statement)


@follows(mkdir('counts_trans.dir'))
@transform(xt_tag,
           regex("mapped_files.dir/(\S+)_XT.bam"),
           r"counts_trans.dir/\1.counts.tsv.gz")
def count_trans(infile, outfile):
    '''Use umi-tools to collapse UMIs and generate counts table'''

    statement = '''umi_tools count --per-gene --gene-tag=XT -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@transform(xt_tag,
           regex("mapped_files.dir/(\S+)_XT.bam"),
           r"counts_trans.dir/\1.counts_unique.tsv.gz")
def count_trans_unique(infile, outfile):
    '''Use umi-tools to count over the umis'''

    statement = '''umi_tools count --per-gene --method=unique --gene-tag=XT -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@merge(count_trans, "counts_trans.dir/counts.tsv.gz")
def merge_count(infiles, outfile):
    '''merge counts from ech sample into one'''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@merge(count_trans_unique, "counts_trans.dir/counts_unique.tsv.gz")
def merge_count_unique(infiles, outfile):
    '''merge counts from ech sample into one'''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


#############################
## Gene level analysis ######
#############################

@transform(tso_umi,
           regex("processed_fastq.dir/(\S+)_polya_tso_UMI.fastq.gz"),
           r"mapped_files.dir/\1_gene.sam")
def mapping_gene(infile, outfile):
    '''map using minimap2 for the geness'''

    if PARAMS['minimap2_splitprefix']:
        statement = """minimap2 -ax splice  -k 14 --split-prefix tmp_sam_ --sam-hit-only --secondary=no --junc-bed %(junc_bed)s %(genome_fasta)s %(infile)s > %(outfile)s  2> %(outfile)s.log"""

    else:
        statement = '''minimap2 -ax splice  -k 14 --sam-hit-only --secondary=no --junc-bed %(junc_bed)s %(genome_fasta)s %(infile)s > %(outfile)s  2> %(outfile)s.log'''

    P.run(statement, job_memory="40G")


@transform(mapping_gene,
           regex("mapped_files.dir/(\S+)_gene.sam"),
           r"mapped_files.dir/\1_gene_sorted.bam")
def samtools_sort(infile, outfile):
    '''strip sequence and then sort the bam file'''

    name = infile.replace("_gene.sam", "")

    statement = '''cgat bam2bam --method=strip-sequence -L strip.log -I %(infile)s -S %(name)s_gene_strip.sam &&
                   samtools view -bh %(name)s_gene_strip.sam > %(name)s_gene.bam &&
                   samtools sort %(name)s_gene.bam -o %(name)s_gene_sorted.bam &&
                   samtools index %(name)s_gene_sorted.bam'''

    P.run(statement)


@transform(samtools_sort,
           regex("mapped_files.dir/(\S+)_gene_sorted.bam"),
           r"mapped_files.dir/\1_featurecounts_gene_sorted.bam")
def featurecounts(infile, outfile):
    '''run featurecounts over the bam file'''

    name = infile.replace("_gene_sorted.bam", "")
    statement = '''featureCounts -a %(gtf)s -o %(name)s_gene_assigned.txt -R BAM %(infile)s &&
                   samtools sort %(infile)s.featureCounts.bam -o %(name)s_featurecounts_gene_sorted.bam &&
                   samtools index %(name)s_featurecounts_gene_sorted.bam'''

    P.run(statement)


@follows(mkdir('counts_genes.dir'))
@transform(featurecounts,
           regex("mapped_files.dir/(\S+)_featurecounts_gene_sorted.bam"),
           r"counts_genes.dir/\1.counts_gene.tsv.gz")
def count_gene(infile, outfile):
    '''Use umi-tools to collapse UMIs and generate counts table'''

    statement = '''umi_tools count --per-gene --gene-tag=XT -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@merge(count_gene, "counts_genes.dir/counts_gene.tsv.gz")
def merge_count_gene(infiles, outfile):
    '''merge counts from ech sample into one'''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@transform(featurecounts,
           regex("mapped_files.dir/(\S+)_featurecounts_gene_sorted.bam"),
           r"counts_genes.dir/\1.count_gene_unique.tsv.gz")
def count_gene_unique(infile, outfile):
    '''Use umi-tools to collapse UMIs and generate counts table'''

    statement = '''umi_tools count --per-gene --gene-tag=XT --method=unique -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@merge(count_gene_unique, "counts_genes.dir/gene_counts_unique.tsv.gz")
def merge_count_gene_unique(infiles, outfile):
    '''merge counts from ech sample into one'''

    df = merge_feature_data(infiles)
    df = df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


########
### analyse without UMI sequences
#########

@follows(featurecounts)
@originate("counts_gene_noumis.tsv.gz")
def merge_featurecounts(outfile):
    ''' '''

    infiles = glob.glob("mapped_files.dir/*_gene_assigned.txt")
    final_df = merge_featurecounts_data(infiles)
    names = [x.replace("_gene_assigned.txt", "") for x in infiles]
    final_df.columns = names
    df = final_df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@follows(merge_count, merge_count_unique, merge_count_gene, merge_count_gene_unique, merge_featurecounts)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
