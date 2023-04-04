"""
===============
Pipeline fusion
===============

Overview
========

This is a Python script that serves as a pipeline to annotate fusion transcripts
in mapped BAM files from a bulk nanopore experiment. The script relies on Ruffus,
a Python library for managing computational pipelines, and the CGAT-core library.


Pipeline tasks
==============

The pipeline performs the following steps:

* Create a BAM file containing only chimeric reads.
* Convert the input BED file to a tabix-indexed file.
* Annotate the fusion genes and the original genes.
* Generate output BED files for the annotated fusion genes.
* Generate counts for each fusion gene.

Input files
===========

The pipeline takes input files in BAM format, mapped using
pipeline_counts.py from a FASTQ file containing nanopore reads
with trimers at the 5' and 3' ends. The input data should be
placed in the data.dir folder.

Pipeline output
===============

The output of the pipeline includes annotated BAM files
and count files for each fusion gene. The script also
generates intermediate files during the process.

Usage
=====

The pipeline requires a configured pipeline.yml file,
which contains various settings and parameters required for
the pipeline to run. 

To generate the config file to change the running of the pipeline you need to
run:

tallytrin fusion config

This will generate a pipeline.yml file that the user can modify to change the
output of the pipeline. Once the user has modified the pipeline.yml file the
pipeline can then be ran using the following commandline command:

tallytrin fusion make full -v5

You can run the pipeline locally (without a cluster) using --local

tallytrin fusion make full -v5 --local


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


SEQUENCESUFFIXES = ("*.bam")

FASTQTARGET = tuple([os.path.join("data.dir/", suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])



@transform("data.dir/*.bam",
         regex("data.dir/(\S+).bam"),
         r"\1_SA.bam")
def make_sabam(infile, outfile):
    '''
    Creates a BAM file containing only chimeric reads
    from the input BAM file.
    '''

    bamfile = pysam.AlignmentFile(infile, "rb")
    split = pysam.AlignmentFile(outfile, "wb", template=bamfile)

    for line in bamfile:
        if line.has_tag("SA"):
            split.write(line)
        else:
            pass

    bamfile.close()
    split.close()


@transform(PARAMS['bed'],
           suffix(".bed"),
           ".bed.gz.tbi")
def tabix_bed(infile, outfile):
    '''
    Converts a BED file into a tabix-indexed file for faster querying.
    '''

    statement = '''cat %(infile)s | bgzip -c > %(infile)s.gz &&
                   tabix -p bed %(infile)s.gz'''

    P.run(statement)



@transform(make_sabam,
           regex("(\S+).bam"),
           add_inputs(tabix_bed),
           r"\1_FusionAnnotate.bam")
def fusion_annotate(infiles, outfile):
    '''
    Annotates fusion and original genes in the input BAM file using the given BED file.
    '''

    infile, bed = infiles
    bed = bed.replace(".tbi", "")

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python  %(PYTHON_ROOT)s/fusion_annotate.py --infile=%(infile)s --outfile=%(outfile)s --bedfile=%(bed)s'''

    P.run(statement, job_options="-t 24:00:00")


@transform(fusion_annotate,
           regex("(\S+)_FusionAnnotate.bam"),
           add_inputs(tabix_bed),
           r"\1_FinalAnnotate.bam")
def gene_annotate(infiles, outfile):
    '''
    Annotates genes and original genes in the input BAM file using the given BED file.
    '''

    infile, bed = infiles
    bed = bed.replace(".tbi", "")
    outfile_tmp = outfile + ".tmp"

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python  %(PYTHON_ROOT)s/gene_annotate.py --infile=%(infile)s --outfile=%(outfile_tmp)s --bedfile=%(bed)s &&
                   samtools sort %(outfile_tmp)s -o %(outfile)s &&
                   rm -rf %(outfile_tmp)s'''

    P.run(statement, job_options="-t 24:00:00")


@transform(gene_annotate,
           regex("(\S+)_FinalAnnotate.bam"),
           [r"\1_fusion1.bed", r"\1_fusion2.bed"])
def generate_bedout(infile, outfiles):
    '''
    Generates output BED files for the annotated fusion genes in the input BAM file.
    '''

    outfile1, outfile2 = outfiles
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/bed_fusion.py --infile=%(infile)s --bed1=%(outfile1)s --bed2=%(outfile2)s'''

    P.run(statement, job_options="-t 24:00:00")


@transform(generate_bedout,
           suffix("_fusion1.bed"),
           "_counts.txt")
def generate_counts(infiles, outfile):
    '''
    Generates counts for each fusion gene using the two input BED files.
    '''

    bed1, bed2 = infiles

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/generate_counts.py --bed1=%(bed1)s --bed2=%(bed2)s --outfile=%(outfile)s'''

    P.run(statement, job_options="-t 24:00:00")


@follows(generate_counts)
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
