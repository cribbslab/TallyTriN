"""
=======================
Pipeline 10X short read
=======================


Overview
==================

This Python code is a pipeline for processing 10X single-cell short read sequencing data using the primer sequence as a common molecular identifier to show the accuracy of single-cell sequencing. The goal of this pipeline is to generate both a gene-level and transcript-level market matrix count matrix (.mtx).

Input files
-----------

The pipeline requires the following inputs:

* paired fastq files named name.fastq.1.gz and name.fastq.2.gz
* a transcriptome genome of your choice
* a genome fasta file
* a gtf file

Pipeline Tasks
==============

The pipeline performs the following tasks:

* 


Pipeline output
===============

The two main outputs of the pipeline are:

* A .mtx file for the gene-level analysis, which can be found within the directory mtx.dir/


These outputs can be easily imported into R using the read_count_output() function within library(BUSpaRse).


Usage
=====

To generate the config file to change the running of the pipeline you need to
run:

tallytrin 10x_shortread config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin 10x_shortread make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin 10x_shortread make full -v5 --local



Code
==================

"""
from ruffus import *

import sys
import os
import re
import sqlite3
import glob

import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.database as database

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
        DATADIR = "data"
    else:
        DATADIR = PARAMS['data']


def connect():
    ''' Connect to database'''

    dbh = sqlite3.connect('csvdb')

    return dbh


SEQUENCESUFFIXES = ("*.fastq.1.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])


@transform(SEQUENCEFILES,
           regex("data.dir/(\S+).fastq.1.gz"),
           r"\1.whitelist.txt")
def identify_barcode(infile, outfile):
    '''
    '''

    barcode_pattern = PARAMS['umi_tools_pattern']
    cells = PARAMS['cells']
    outfile = outfile.replace('data.dir/', '')
    statement = '''umi_tools whitelist --stdin %(infile)s
                                    --bc-pattern=%(barcode_pattern)s
                                    --set-cell-number=%(cells)s
                                    --log2stderr > %(outfile)s'''

    P.run(statement, job_options='-t 24:00:00')


@follows(mkdir('extracted_fastq.dir'))
@transform(SEQUENCEFILES,
           regex('data.dir/(\S+).fastq.1.gz'),
           add_inputs(identify_barcode),
           r'extracted_fastq.dir/\1_extracted.fastq.1.gz')
def extract_bcumi(infiles, outfile):

    infile, whitelist = infiles
    read2_in = infile.replace('.fastq.1.gz', '.fastq.2.gz')
    read2 = outfile.replace('.fastq.1.gz', '.fastq.2.gz')
    barcode_pattern = PARAMS['umi_tools_pattern']    

    statement = '''umi_tools extract --bc-pattern=%(barcode_pattern)s
                   --stdin %(infile)s
                   --stdout %(outfile)s
                   --read2-in %(read2_in)s
                   --read2-out %(read2)s
                   --whitelist=%(whitelist)s'''

    P.run(statement, job_options='-t 24:00:00')


@follows(mkdir("hisat_mapping.dir"))
@transform(extract_bcumi,
           regex('extracted_fastq.dir/(\S+)_extracted.fastq.1.gz'),
           r'hisat_mapping.dir/\1.bam')
def hisat2_map(infile, outfile):
    '''
    '''

    read1 = "".join(infile)
    read2 = read1.replace('.1.gz', '.2.gz')

    job_threads = PARAMS['hisat2_threads']
    index_prefix = PARAMS['hisat2_index_dir'] + '/' + PARAMS['hisat2_index_name']

    statement = '''ln -s %(read2)s tmp.fastq.gz && hisat2 %(hisat2_options)s -x %(index_prefix)s --threads %(job_threads)i
                   -U tmp.fastq.gz 2> %(outfile)s.log | samtools view -bS > %(outfile)s.tmp
                   2>> %(outfile)s.log && samtools sort -o %(outfile)s %(outfile)s.tmp &&
                   rm -rf %(outfile)s.tmp'''

    P.run(statement, job_threads=job_threads, job_options='-t 24:00:00',
          job_memory=PARAMS['hisat2_memory'])


@transform(hisat2_map,
           regex('hisat_mapping.dir/(\S+).bam'),
           r'\1.bam.featureCounts.bam')
def featurecounts(infile, outfile):
    '''
    '''

    feature_threads = PARAMS['featurecounts_threads']
    geneset = PARAMS['featurecounts_gtf']

    statement = '''featureCounts -a %(geneset)s -o gene_assigned 
                   -R BAM %(infile)s -T %(feature_threads)s '''

    P.run(statement, job_threads=feature_threads)


@transform(featurecounts,
           regex('(\S+).featureCounts.bam'),
           r'\1_sorted.bam')
def sort_index(infile, outfile):
    '''
    '''

    statement = '''samtools sort %(infile)s -o %(outfile)s &&
                   samtools index %(outfile)s'''

    P.run(statement)


@transform(sort_index,
           regex("(\S+).bam_sorted.bam"),
           r"counts.tsv.gz")
def count(infile, outfile):
    '''Counts the reads in the input BAM file using umi_tools with
    unique method, per gene and per cell. The output is a compressed
    TSV file called "counts.tsv.gz".'''

    statement = '''umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell  -I %(infile)s -S counts.tsv.gz'''

    P.run(statement, job_options='-t 24:00:00')


@follows(mkdir("mtx.dir"))
@transform(count,
           regex("counts.tsv.gz"),
           r"mtx.dir/genes.mtx")
def convert_tomtx(infile, outfile):
    '''
    Converts the count matrix in the input TSV file to a .mtx format.
    The output is stored in the "mtx.dir" directory.
    '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/save_mtx.py --data=%(infile)s --dir=mtx.dir/'''

    P.run(statement, job_memory="100G", job_options='-t 24:00:00')


@follows(convert_tomtx)
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
