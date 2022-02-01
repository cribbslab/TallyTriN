"""===========================
Pipeline count
===========================

Overview
========

The aim of this pipeline is to take a mapped bam file from a bulk nanopore experiment and then annotate the fusion transcripts.

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

.bam file mapped using pipeline_counts.py using a fastq.gz file of nanopore reads that have been sequenced with trimers at
the 5' and 3' end. Data should be added to the data.dir folder.

Pipeline output
===============



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
    '''Create a bam file with all of the reads listed as chimeric'''

    E.warn(infile)
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
    ''' convert bed file to tabix'''

    statement = '''cat %(infile)s | bgzip -c > %(infile)s.gz &&
                   tabix -p bed %(infile)s.gz'''

    P.run(statement)



@transform(make_sabam,
           regex("(\S+).bam"),
           add_inputs(tabix_bed),
           r"\1_FusionAnnotate.bam")
def fusion_annotate(infiles, outfile):
    '''Annotate fusion and original gene'''

    infile, bed = infiles
    bed = bed.replace(".tbi", "")

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python  %(PYTHON_ROOT)s/fusion_annotate.py --infile=%(infile)s --outfile=%(outfile)s --bedfile=%(bed)s'''

    P.run(statement)


@transform(fusion_annotate,
           regex("(\S+)_FusionAnnotate.bam"),
           add_inputs(tabix_bed),
           r"\1_FinalAnnotate.bam")
def gene_annotate(infiles, outfile):
    '''Annotate gene and original gene'''

    infile, bed = infiles
    bed = bed.replace(".tbi", "")
    outfile_tmp = outfile + ".tmp"

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python  %(PYTHON_ROOT)s/gene_annotate.py --infile=%(infile)s --outfile=%(outfile_tmp)s --bedfile=%(bed)s &&
                   samtools sort %(outfile_tmp)s -o %(outfile)s &&
                   rm -rf %(outfile_tmp)s'''

    P.run(statement)


@transform(gene_annotate,
           regex("(\S+)_FinalAnnotate.bam"),
           [r"\1_fusion1.bed", r"\1_fusion2.bed"])
def generate_bedout(infile, outfiles):
    '''Generate output bed file '''

    outfile1, outfile2 = outfiles
    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/bed_fusion.py --infile=%(infile)s --bed1=%(outfile1)s --bed2=%(outfile2)s'''

    P.run(statement)


@transform(generate_bedout,
           suffix("_fusion1.bed"),
           "_counts.txt")
def generate_counts(infiles, outfile):
    '''Generate counts for each fusion gene'''

    bed1, bed2 = infiles

    PYTHON_ROOT = os.path.join(os.path.dirname(__file__), "python/")

    statement = '''python %(PYTHON_ROOT)s/generate_counts.py --bed1=%(bed1)s --bed2=%(bed2)s --outfile=%(outfile)s'''

    P.run(statement)


@transform(generate_counts,
           suffix("_counts.txt"),
           "_finalcounts.txt")
def generate_finalcounts(infile, outfile):
    '''Generate final counts for each fusion gene irrespective of the orientation of the fusion'''

    R_ROOT = os.path.join(os.path.dirname(__file__), "R/")

    statement = '''Rscript %(R_ROOT)s/merge_counts.R --input=%(infile)s --out=%(outfile)s'''

    P.run(statement)


@follows(generate_finalcounts)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
