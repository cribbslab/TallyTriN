"""
=================
Pipeline illumina
=================

Overview
==================

This pipeline processes Illumina sequenced FASTQ files with trimer UMI
(Unique Molecular Identifiers) for single-cell analysis using
hisat2 and featurecounts. The pipeline processes the input FASTQ files,
corrects the double barcodes/UMI in read1, and generates a collapsed
FASTQ file compatible with downstream workflows.

Pipeline tasks
==============

The code consists of several functions, each handling a specific step in the processing and analysis of the input FASTQ files.

* UMI Correction: The correct_umis() function corrects the UMI sequences in the input Illumina data by processing read1 and read2 FASTQ files.
* Read Mapping: The map_hisat2() function maps the reads to a reference genome using HISAT2, a fast and efficient read aligner.
* SAM to BAM Conversion and Sorting: The run_samtools() function converts the output SAM files from the read mapping step into BAM format, sorts them, and indexes the sorted BAM files.
* Feature Counting: The featurecounts() function runs featureCounts to count the number of reads associated with each genomic feature (e.g., gene) and outputs a BAM file.
* Gene Counting with UMIs: The count_genes() function uses UMI-tools to count the reads for each gene, considering the UMIs.
* Gene Counting without UMIs: The count_genes_noumis() function uses UMI-tools to count the reads for each gene without considering the UMIs (unique method).
* Merging Feature Count Data: The merge_featurecounts_data() function merges the input files from featureCounts and returns a DataFrame with the counts for each gene.
* Merging Gene Count Data with UMIs: The merge_genes() function merges the gene count data from featureCounts with UMI consideration into a single output file.
* Merging Gene Count Data without UMIs: The merge_genes_noumi() function merges the gene count data from featureCounts without UMI consideration into a single output file.
* Merging Feature Count Data without UMIs: The merge_featurecounts() function merges the output files from featureCounts and generates a single output file with gene counts without UMI consideration.
* Completion: The full() function represents the completion of the pipeline.

Input files
===========

The workflow requires the following input:
* paired gz fastq files following Illumina sequencing

Pipeline output
==================

The output of this pipeline is a "perfect_collapsed.fastq.[1-2].gz" fastq file.
this file has been barcode corrected an collapsed into single nucleotides so that
it is compatible with downstream workflows.    

Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallytrin illumina config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin illumina make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin illumina make full -v5 --local


Code
==================

"""
from ruffus import *

import sys
import os
import re
import sqlite3
import glob
import pandas as pd

import cgatcore.pipeline as P
import cgatcore.experiment as E

# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


# Determine the location of the input fastq files
try:
    PARAMS['data']
except NameError:
    DATADIR = "."
else:
    if PARAMS['data'] == 0:
        DATADIR = "."
    elif PARAMS['data'] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS['data']


SEQUENCESUFFIXES = ("*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])



@follows(mkdir("corrected_umis.dir"))
@transform(SEQUENCEFILES,
           regex("data.dir/(\S+).fastq.1.gz"),
           r"corrected_umis.dir/\1_corrected.fastq.gz")
def correct_umis(infile, outfile):
    '''
    Corrects the UMI sequences in the Illumina data by
    processing read1 and read2 FASTQ files.
    '''

    read1 = infile
    read2 = infile.replace(".fastq.1.gz", ".fastq.2.gz")

    name = infile.replace("", ".fastq.1.gz")
    name = name.replace("data.dir/", "")

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")


    if PARAMS["correct"]:
        statement = '''python %(PYTHON_ROOT)s/correct_illumina_umi.py --read1=%(read1)s --read2=%(read2)s --outname=%(outfile)s --errors=%(errors)s'''
    else:
        statement = '''python %(PYTHON_ROOT)s/uncorrect_illumina.py --read1=%(read1)s --read2=%(read2)s --outname=%(outfile)s'''


    P.run(statement)


@transform(correct_umis,
           regex("corrected_umis.dir/(\S+)_corrected.fastq.gz"),
           r"mapped.dir/\1.sam")
def map_hisat2(infile, outfile):
    """
    Maps the reads using HISAT2.
    """

    statement = """hisat2 -x %(hisat2_index)s -U %(infile)s -S %(outfile)s"""

    P.run(statement)


@transform(map_hisat2,
           regex("mapped.dir/(\S+).sam"),
           r"mapped.dir/\1_sorted.bam")
def run_samtools(infile, outfile):
    """
    Converts SAM files to BAM, sorts them, and indexes the sorted BAM files.
    """

    statement = """samtools view -bh  %(infile)s > %(infile)s_final_gene.bam &&
                   samtools sort %(infile)s_final_gene.bam -o %(outfile)s &&
                   samtools index %(outfile)s"""
    
    P.run(statement)


@transform(run_samtools,
           regex("mapped.dir/(\S+)_sorted.bam"),
           r"mapped.dir/\1_Aligned_final_gene_sorted.bam")
def featurecounts(infile, outfile):
    """
    Runs featureCounts to count the number of reads associated
    with each genomic feature (e.g., gene) and outputs a BAM file.
    """

    name = outfile.replace("_Aligned_final_gene_sorted.bam", "_gene_assigned")

    statement = """featureCounts -a %(gtf)s -o %(name)s -R BAM %(infile)s &&
                   samtools sort %(infile)s.featureCounts.bam  -o %(outfile)s &&
                   samtools index %(outfile)s"""

    P.run(statement)


@follows(mkdir("featurecounts.dir"))
@transform(featurecounts,
           regex("mapped.dir/(\S+)_Aligned_final_gene_sorted.bam"),
           r"featurecounts.dir/\1_counts_genes.tsv.gz")
def count_genes(infile, outfile):
    '''
    Uses UMI-tools to count the reads for each gene, considering the UMIs.
    '''

    statement = '''umi_tools count --per-gene --gene-tag=XT  -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@transform(featurecounts,
           regex("mapped.dir/(\S+)_Aligned_final_gene_sorted.bam"),
           r"featurecounts.dir/\1_counts_genes_noumis.tsv.gz")
def count_genes_noumis(infile, outfile):
    '''
    Uses UMI-tools to count the reads for each gene without considering the UMIs (unique method).
    '''

    statement = '''umi_tools count --method=unique --per-gene --gene-tag=XT  -I %(infile)s -S %(outfile)s'''

    P.run(statement)


def merge_featurecounts_data(infiles):
    '''
    Merges the input files from featureCounts and returns a DataFrame with the counts for each gene.
    '''

    final_df = pd.DataFrame()
    for infile in infiles:
    
        tmp_df = pd.read_table(infile, sep="\t", header=0, index_col=0, skiprows = 0, compression='gzip')
        tmp_df = tmp_df.iloc[:,-1:]
        tmp_df.columns = ["count"]
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True, suffixes=("","_drop"))

    names = [x.replace("", "") for x in infiles]
    final_df.columns = names
    return final_df

@follows(count_genes)
@originate("counts_gene.tsv.gz")
def merge_genes(outfile):
    '''
    Merges the gene count data from featureCounts with UMI consideration into a single output file.
    '''

    infiles = glob.glob("featurecounts.dir/*_counts_genes.tsv.gz")
    final_df = merge_featurecounts_data(infiles)
    names = [x.replace("_counts_genes.tsv.gz", "") for x in infiles]
    final_df.columns = names
    df = final_df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@follows(count_genes_noumis)
@originate("counts_gene_unique.tsv.gz")
def merge_genes_noumi(outfile):
    '''
    Merges the gene count data from featureCounts without UMI consideration into a single output file.
    '''

    infiles = glob.glob("featurecounts.dir/*_counts_genes_noumis.tsv.gz")
    final_df = merge_featurecounts_data(infiles)
    names = [x.replace("_counts_genes_noumis.tsv.gz", "") for x in infiles]
    final_df.columns = names
    df = final_df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")


@follows(featurecounts)
@originate("counts_gene_noumis.tsv.gz")
def merge_featurecounts(outfile):
    '''
    Merges the output files from featureCounts and generates a single output
    file with gene counts without UMI consideration.
    '''

    infiles = glob.glob("mapped.dir/*_gene_assigned")
    final_df = pd.DataFrame()

    for infile in infiles:
    
        tmp_df = pd.read_table(infile, sep="\t", header=0, index_col=0, skiprows = 1)
        tmp_df = tmp_df.iloc[:,-1:]
        tmp_df.columns = ["count"]
        final_df = final_df.merge(tmp_df, how="outer", left_index=True, right_index=True, suffixes=("","_drop"))
    names = [x.replace("_gene_assigned", "") for x in infiles]
    final_df.columns = names
    
    df = final_df.fillna(0)
    df.to_csv(outfile, sep="\t", compression="gzip")

@follows(merge_genes, merge_genes_noumi, merge_featurecounts)
def full():
    '''
    A placeholder function that serves as a checkpoint
    to run all previous ruffus tasks and ensure that all
    previous tasks are completed.
    '''
    pass


def main(argv=None):
    '''
    The main function that runs the pipeline using the cgatcore.pipeline module.
    Takes an optional argument list (default is sys.argv).

    Please note that some of these functions use external Python scripts or
    tools. For a complete understanding of their functionality, it is
    necessary to examine the code of those scripts as well.
    '''
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
