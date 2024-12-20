"""
============================
Pipeline single-cell_hydrop
============================

Overview
==================

This code is for a Python script that processes single-cell sequencing
data from long-read monomer barcode macosko sequenced using nanopore. 
It generates a gene-level and transcript-level counts matrix in
market matrix format.

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

A single fastq file generated by guppy basecalling.
A transcriptome genome of your choice.
A genome fasta file.
A minimap2 junction bed generated following minimap2 instructions.
A GTF file.
The main outputs of the pipeline include:

Output files
============

A market matrix format (.mtx) output of the counts for transcripts, found within the mtx.dir/ directory.
A .mtx file for the gene-level analysis, found within the mtx_gene.dir/ directory.

These outputs can be imported into R using the read_count_output() function within the library(BUSpaRse) package.

Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallytrin singlecell_macosko config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin singlecell_macosko make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin singlecell_macosko make full -v5 --local


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



SEQUENCESUFFIXES = ("*.fastq")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])


@follows(mkdir("split_tmp.dir"))
@split(DATADIR + '/*.fastq.gz', "split_tmp.dir/out*")
def split_fastq(infile, outfiles):
    '''
    Splits the input FASTQ file into smaller pieces using the seqkit split
    command. This function takes an input FASTQ file and creates multiple
    smaller FASTQ files based on the desired split size specified in the
    configuration.
    '''

    infile = "".join(infile)

    statement = '''zcat %(infile)s | split -l %(split)s - out. &&
                   mv out*.* split_tmp.dir/'''

    P.run(statement)


@follows(mkdir("polyA_correct.dir"))
@transform(split_fastq,
           regex("split_tmp.dir/out.(\S+)"),
           r"polyA_correct.dir/\1_correct_polya.fastq")
def correct_polyA(infile, outfile):
    '''
    Performs quality control on the split FASTQ files by correcting
    the polyA tail and identifying barcode and UMI sequences using
    a custom script correct_polyA.py. This function takes an input
    split FASTQ file and generates a corrected FASTQ file.
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/complement_polyA_singlecell.py --infile=%(infile)s --outname=%(outfile)s'''

    P.run(statement)


@follows(mkdir("polyA_umi.dir"))
@transform(correct_polyA,
           regex("polyA_correct.dir/(\S+)_correct_polya.fastq"),
           r"polyA_umi.dir/\1.fastq.gz")
def identify_bcumi(infile, outfile):
    '''
    Identifies the barcode and unique molecular identifier (UMI) sequences from
    the corrected FASTQ file using a custom script identify_bcumi.py. This
    function takes an input corrected FASTQ file and generates a FASTQ
    file with annotated barcode and UMI information in the read names.
    '''

    name = outfile.replace("polyA_umi.dir/", "")
    name = name.replace(".fastq.gz", "")
    barcode_len = PARAMS['barcode_len']
    cmimode = PARAMS['cmimode']

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/identify_bcumi_hydrop.py --outfile=%(outfile)s --infile=%(infile)s --whitelist=polyA_umi.dir/%(name)s.whitelist.txt --barcode_len=%(barcode_len)s --cmimode=%(cmimode)s'''

    P.run(statement)



@merge(identify_bcumi, "whitelist.txt")
def merge_whitelist(infiles, outfile):
    '''
    merge whitelists
    '''

    whitelists = []

    for i in infiles:

        whitelists.append(i.replace(".fastq.gz", ".whitelist.txt"))

    whitelist_files = " ".join(whitelists)

    statement = '''cat %(whitelist_files)s | sort | uniq > %(outfile)s'''

    P.run(statement)


@follows(merge_whitelist)
@follows(mkdir("correct_reads.dir"))
@transform(identify_bcumi,
           regex("polyA_umi.dir/(\S+).fastq.gz"),
           r"correct_reads.dir/\1.fastq.gz")
def correct_reads(infile, outfile):
    '''Correct the barcodes using majority vote'''

    infile = infile
    cells = PARAMS['cells']

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")


    statement = '''python %(PYTHON_ROOT)s/correct_10xbarcode.py --infile=%(infile)s --outfile=%(outfile)s --cells=%(cells)s --whitelist=whitelist.txt  --cmimode=0 --umi=8'''

    P.run(statement, job_options='-t 24:00:00')


@merge(correct_reads, "merge_corrected.fastq.gz")
def merge_correct_reads(infiles, outfile):
    '''Merge the corrected reads '''

    infile = []

    for i in infiles:
        infile.append(str(i))

    infile_join = " ".join(infile)



    statement = '''cat %(infile_join)s > %(outfile)s'''

    P.run(statement)


@transform(merge_correct_reads,
           regex("merge_corrected.fastq.gz"),
           r"final.sam")
def mapping(infile, outfile):
    '''
    Maps the reads from the merged FASTQ file to a transcriptome genome using
    the minimap2 aligner. This function takes an input merged FASTQ file
    and generates a SAM (Sequence Alignment/Map) file containing the
    read alignments.
    '''


    cdna = PARAMS['minimap2_fasta_cdna']
    options = PARAMS['minimap2_options']

    statement = '''minimap2  %(options)s %(cdna)s  %(infile)s > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement, job_options='-t 24:00:00')


@transform(mapping,
           regex("final.sam"),
           r"final_sorted.bam")
def run_samtools(infile, outfile):
    '''convert sam to bam and sort -F 272'''

    statement = '''samtools view -bS %(infile)s > final.bam &&
                   samtools sort final.bam -o final_sorted.bam &&
                   samtools index final_sorted.bam'''

    P.run(statement, job_options='-t 24:00:00')


@transform(run_samtools,
           regex("final_sorted.bam"),
           r"final_XT.bam")
def add_xt_tag(infile, outfile):
    '''
    Adds the transcript name to the XT tag in the BAM file using a custom
    script add_xt_tag.py. This function takes an input sorted BAM file
    and generates a BAM file with the added XT tag containing
    transcript information.
    '''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/xt_tag_nano.py --infile=%(infile)s --outfile=%(outfile)s &&
                   samtools index %(outfile)s'''

    P.run(statement)


@transform(add_xt_tag,
           regex("final_XT.bam"),
           r"counts.tsv.gz")
def count(infile, outfile):
    '''
    Counts the reads in the BAM file with XT tags using the umi_tools package.
    The counting is performed per gene and using a unique method.
    '''

    statement = '''umi_tools count --per-gene --gene-tag=XT --per-cell  -I %(infile)s -S counts.tsv.gz'''

    P.run(statement)


@follows(mkdir("mtx.dir"))
@transform(count,
           regex("counts.tsv.gz"),
           r"mtx.dir/genes.mtx")
def convert_tomtx(infile, outfile):
    '''
    This function takes an input BAM file with XT tags and generates a counts matrix file in the MTX format.
    '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/save_mtx.py --data=%(infile)s --dir=mtx.dir/'''

    P.run(statement, job_memory="250G")


@merge(identify_bcumi, "merge_uncorrected.fastq.gz")
def merge_uncorrect_reads(infiles, outfile):
    '''Merge the reads that are still containing trimer barcode and umi'''

    infile = []

    for i in infiles:
        infile.append(str(i))

    infile_join = " ".join(infile)



    statement = '''cat %(infile_join)s > %(outfile)s'''

    P.run(statement)


@transform(merge_uncorrect_reads,
           regex("merge_uncorrected.fastq.gz"),
           r"final_uncorrected.sam")
def mapping_trimer(infile, outfile):
    '''Run minimap2 to map the fastq file'''


    cdna = PARAMS['minimap2_fasta_cdna']
    options = PARAMS['minimap2_options']

    statement = '''minimap2  %(options)s %(cdna)s  %(infile)s > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement, job_options='-t 24:00:00')


@follows(mkdir("collapse_reads.dir"))
@transform(identify_bcumi,
           regex("polyA_umi.dir/(\S+).fastq.gz"),
           r"collapse_reads.dir/\1.fastq.gz")
def collapse_reads(infile, outfile):
    '''Correct the barcodes by picking first ucleotide in the barcode and umi'''

    infile = infile

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    if PARAMS['trimer_beads']:

        statement = '''python %(PYTHON_ROOT)s/single_nucleotide_select.py --infile=%(infile)s --outfile=%(outfile)s'''

    else:
        statement = '''cp %(infile)s %(outfile)s '''


    P.run(statement)


@merge(collapse_reads, "merge_collapsed.fastq.gz")
def merge_singlenuc_reads(infiles, outfile):
    '''Merge the reads that were collapsed into single nucleotides'''

    infile = []

    for i in infiles:
        infile.append(str(i))

    infile_join = " ".join(infile)



    statement = '''cat %(infile_join)s > %(outfile)s'''

    P.run(statement)


@transform(merge_singlenuc_reads,
           regex("merge_collapsed.fastq.gz"),
           r"final_collapsed.sam")
def mapping_collapsed(infile, outfile):
    '''Run minimap2 to map the fastq file'''


    cdna = PARAMS['minimap2_fasta_cdna']
    options = PARAMS['minimap2_options']

    statement = '''minimap2  %(options)s %(cdna)s  %(infile)s > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement, job_options='-t 24:00:00')



@transform(mapping_collapsed,
           regex("final_collapsed.sam"),
           r"final_sorted_collapsed.bam")
def run_samtools_collapsed(infile, outfile):
    '''convert sam to bam and sort -F 272'''

    statement = '''samtools view -bS %(infile)s > final_collapsed.bam &&
                   samtools sort final_collapsed.bam -o final_sorted_collapsed.bam &&
                   samtools index final_sorted_collapsed.bam'''

    P.run(statement, job_options='-t 24:00:00')


@transform(run_samtools_collapsed,
           regex("final_sorted_collapsed.bam"),
           r"final_XT_collapsed.bam")
def add_xt_tag_collapsed(infile, outfile):
    '''Add trancript name to XT tag in bam file so umi-tools counts can be  perfromed'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/xt_tag_nano.py --infile=%(infile)s --outfile=%(outfile)s &&
                   samtools index %(outfile)s'''

    P.run(statement)


@transform(add_xt_tag_collapsed,
           regex("final_XT_collapsed.bam"),
           r"counts_collapsed.tsv.gz")
def count_collapsed(infile, outfile):
    '''use umi_tools to count the reads - need to adapt umi tools to double oligo'''

    statement = '''umi_tools count --method unique --per-gene --gene-tag=XT --per-cell  -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@follows(mkdir("mtx_collapsed.dir"))
@transform(count_collapsed,
           regex("counts_collapsed.tsv.gz"),
           r"mtx_collapsed.dir/genes.mtx")
def convert_tomtx_collapsed(infile, outfile):
    ''' '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/save_mtx.py --data=%(infile)s --dir=mtx_collapsed.dir/ --filter=3'''

    P.run(statement, job_memory="250G")



@transform(add_xt_tag_collapsed,
           regex("final_XT_collapsed.bam"),
           r"counts_collapsed_directional.tsv.gz")
def count_collapsed_direction(infile, outfile):
    '''use umi_tools to count the reads - need to adapt umi tools to double oligo'''

    statement = '''umi_tools count  --per-gene --gene-tag=XT --per-cell  -I %(infile)s -S %(outfile)s'''

    P.run(statement)


@follows(mkdir("mtx_collapsed_directional.dir"))
@transform(count_collapsed_direction,
           regex("counts_collapsed_directional.tsv.gz"),
           r"mtx_collapsed_directional.dir/genes.mtx")
def convert_tomtx_directional(infile, outfile):
    ''' '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/save_mtx.py --data=%(infile)s --dir=mtx_collapsed_directional.dir/ --filter=3'''

    P.run(statement, job_memory="250G")


###########################################################################
# Correct the UMIs using greedy
###########################################################################

# Need to input the bam file without collapsing trimers - minimap with fastq before collapsing then
# running ILP


@merge(correct_reads, "merge_trimer.fastq.gz")
def merge_trimer_bcumi(infiles, outfile):
    '''Merge the fastq reads with uncollapsed trime reads '''

    infile = []

    for i in infiles:
        infile.append(str(i))

    infile_join = " ".join(infile)



    statement = '''cat %(infile_join)s > %(outfile)s'''

    P.run(statement)


@transform(merge_trimer_bcumi,
           regex("merge_trimer.fastq.gz"),
           r"final_trimer.sam")
def run_minimap2_trimer(infile, outfile):
    '''Run minimap2 using fastq files with trimer UMIs'''  

    cdna = PARAMS['minimap2_fasta_cdna']
    options = PARAMS['minimap2_options']

    statement = '''minimap2  %(options)s %(cdna)s  %(infile)s > %(outfile)s 2> %(outfile)s.log'''

    P.run(statement, job_options='-t 24:00:00')


@transform(run_minimap2_trimer,
           regex("final_trimer.sam"),
           r"final_sorted_trimer.bam")
def run_samtools_trimer(infile, outfile):
    '''convert sam to bam and sort -F 272'''

    statement = '''samtools view -bS %(infile)s > final_trimer.bam &&
                   samtools sort final_trimer.bam -o final_sorted_trimer.bam &&
                   samtools index final_sorted_trimer.bam'''

    P.run(statement, job_options='-t 24:00:00')


@transform(run_samtools_trimer,
           regex("final_sorted_trimer.bam"),
           r"final_XT_trimer.bam")
def add_xt_tag_trimer(infile, outfile):
    '''Add trancript name to XT tag in bam file so umi-tools counts can be  perfromed'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/xt_tag_nano.py --infile=%(infile)s --outfile=%(outfile)s &&
                   samtools index %(outfile)s'''

    P.run(statement)



@transform(add_xt_tag_trimer,
           regex("final_XT_trimer.bam"),
           r"greedy.csv")
def run_greedy(infile, outfile):
    '''Run greedy algorithm to collapse the UMIs'''

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/greedy_sc.py count -i %(infile)s -t XT -o %(outfile)s'''

    P.run(statement)


@follows(mkdir("mtx_greedy.dir"))
@transform(run_greedy,
           regex("greedy.csv"),
           r"mtx_greedy.dir/genes.mtx")
def convert_tomtx_greedy(infile, outfile):
    ''' '''
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/save_mtx.py --data=%(infile)s --dir=mtx_greedy.dir/'''

    P.run(statement, job_memory="250G")


@follows(convert_tomtx_directional, convert_tomtx_collapsed, convert_tomtx, convert_tomtx_greedy)
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




