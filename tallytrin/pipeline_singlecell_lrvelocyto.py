"""
============================
Pipeline single-cell_macosko
============================

Overview
==================

This code is for a Python script that processes single-cell sequencing
data and generated veloyto output for running downstream analysis at the gene and 
transcript level. It requires you to run the singlecell pipeline and 
generate the mapping_genome.dir/ sorted bam file.

The pipeline uses the CGAT-core library for pipeline management and the
Ruffus library for task management.


Input steps
===========

The pipeline performs the following steps:

Reads the pipeline configuration from the pipeline.yml file.
Splits the input fastq file into smaller pieces.


Input files
===========

To run the pipeline, you will need the following input files:

A genome mapped sorted bam file from running the singlecell pipleine
A minimap2 junction bed generated following minimap2 instructions.
A BED file.

The main outputs of the pipeline include:

Output files
============

A loom file for import into tools such as scVelo

Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallytrin singlecell_levelocyto config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin singlecell_levelocyto make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin singlecell_levelocyto make full -v5 --local


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



SEQUENCESUFFIXES = ("*.bam")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])


@follows(mkdir("unspliced.dir"))
@transform(DATADIR + '/*.bam',
           regex(DATADIR + '/(\S+).bam'),
           r"unspliced.dir/intron_overlapping.bam")
def split_unspliced(infile, outfile):
    '''
    Splits the input BAM file into reads that overlap a intron so
    data can be parsed faster and easier.
    '''

    infile = "".join(infile)

    statement = '''
                bedtools intersect -a %(infile)s -b %(intron_bed)s -f 0.8 -wa -u > %(outfile)s'''

    P.run(statement)


@follows(mkdir("spliced.dir"))
@transform(DATADIR + '/*.bam',
           regex(DATADIR + '/(\S+).bam'),
           r"spliced.dir/intron_nonoverlapping.bam")
def split_spliced(infile, outfile):
    '''
    Splits the input BAM file into reads that overlap a intron so
    data can be parsed faster and easier.
    '''

    infile = "".join(infile)

    statement = '''
                bedtools intersect -a %(infile)s -b %(intron_bed)s -f 0.8 -wa -v > %(outfile)s'''

    P.run(statement)


@transform(split_spliced,
           regex('spliced.dir/(\S+).bam'),
           r"spliced.dir/\1.fastq.gz")
def spliced_fastq(infile, outfile):
    '''
    '''

    infile = "".join(infile)

    statement = '''
                bedtools bamtofastq -i %(infile)s -fq /dev/stdout | gzip > %(outfile)s
                '''

    P.run(statement)


@transform(split_unspliced,
           regex('unspliced.dir/(\S+).bam'),
           r"unspliced.dir/\1.fastq.gz")
def unspliced_fastq(infile, outfile):
    '''
    '''

    infile = "".join(infile)

    statement = '''
                bedtools bamtofastq -i %(infile)s -fq /dev/stdout | gzip > %(outfile)s
                '''

    P.run(statement)


@transform(spliced_fastq,
           regex("spliced.dir/(\S+).fastq.gz"),
           r"spliced.dir/\1.sam")
def mapping_spliced(infile, outfile):
    '''
    Maps the reads from the spliced FASTQ file to a transcriptome genome using
    the minimap2 aligner. This function takes an input merged FASTQ file
    and generates a SAM (Sequence Alignment/Map) file containing the
    read alignments.
    '''


    cdna = PARAMS['minimap2_fasta_cdna']
    options = PARAMS['minimap2_options']

    statement = '''minimap2  %(options)s %(cdna)s  %(infile)s > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement, job_options='-t 24:00:00')



@transform(unspliced_fastq,
           regex("unspliced.dir/(\S+).fastq.gz"),
           r"unspliced.dir/\1.sam")
def mapping_unspliced(infile, outfile):
    '''
    Maps the reads from the unspliced FASTQ file to a transcriptome genome using
    the minimap2 aligner. This function takes an input merged FASTQ file
    and generates a SAM (Sequence Alignment/Map) file containing the
    read alignments.
    '''


    cdna = PARAMS['minimap2_fasta_cdna']
    options = PARAMS['minimap2_options']

    statement = '''minimap2  %(options)s %(cdna)s  %(infile)s > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement, job_options='-t 24:00:00')


@transform(mapping_spliced,
           regex("(\S+).sam"),
           r"spliced.dir/final_trans_sorted.bam")
def run_samtools_spliced_trans(infile, outfile):
    '''convert sam to bam and sort -F 272'''

    statement = '''samtools view -bS %(infile)s > spliced.dir/final_trans.bam &&
                   samtools sort spliced.dir/final_trans.bam -o spliced.dir/final_trans_sorted.bam &&
                   samtools index spliced.dir/final_trans_sorted.bam'''

    P.run(statement, job_options='-t 24:00:00')


@transform(mapping_unspliced,
           regex("(\S+).sam"),
           r"unspliced.dir/final_trans_sorted.bam")
def run_samtools_unspliced_trans(infile, outfile):
    '''convert sam to bam and sort -F 272'''

    statement = '''samtools view -bS %(infile)s > unspliced.dir/final_trans.bam &&
                   samtools sort unspliced.dir/final_trans.bam -o unspliced.dir/final_trans_sorted.bam &&
                   samtools index unspliced.dir/final_trans_sorted.bam'''

    P.run(statement, job_options='-t 24:00:00')


@follows(run_samtools_spliced_trans, run_samtools_unspliced_trans)
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




