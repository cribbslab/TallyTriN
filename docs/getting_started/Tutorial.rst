.. _getting_started-Tutorial:


=============================
Running a pipeline - Tutorial
=============================


Before beginning this tutorial make sure you have the TallyTrin installed correctly,
please see here (see :ref:`getting_started-Installation`) for installation instructions.

As a tutorial example of how to run a TallyTrin workflow we will run the count pipeline.

This workflow is for generating a count matrix for
downstream differential expression analysis using nanopore reads.
The pipeline takes an input fastq file, processes it, and outputs 
a count matrix with samples as columns and rows as either transcripts
or genes. The pipeline makes use of multiple Python libraries and tools
like Minimap2, Samtools, UMI-tools, and mclumi.


Tutorial start
--------------


**1.** First download the tutorial data::

   mkdir count
   cd count
   wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
   wget https://datashare.molbiol.ox.ac.uk/public/project/cribbslab/acribbs/tallytrin_tutorial/hg38.fasta
   wget https://datashare.molbiol.ox.ac.uk/public/project/cribbslab/acribbs/tallytrin_tutorial/geneset_all.bed
   wget https://datashare.molbiol.ox.ac.uk/public/project/cribbslab/acribbs/tallytrin_tutorial/geneset_all.gtf.gz

   gunzip geneset_all.gtf.gz
   gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

You will need to download the test data file from GEO using sra-toolkit (https://github.com/ncbi/sra-tools) here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX18398120&o=acc_s%3Aa


**2.** Next we will generate a configuration yml file so the pipeline output can be modified::

   conda activate tallytrin

   # To show all available pipelines:
   tallytrin -h

   # Generate config
   tallytrin count config

**3.** Modify the config if required::
 
At this stage you would normally modify the config, but in this case the defaults should be fine in 
this case::

  # Config file for pipeline_count.py

  ## general options

  # Copyright statement
  copyright: cribbslab (2021)

  cdna_fasta: Homo_sapiens.GRCh38.cdna.all.fa

  genome_fasta: hg38.fasta
 
  junc_bed: hg38_geneset_all.bed

  gtf: hg38_geneset_all.gtf

  # Specify if the pipeline should run umi correction or not
  correct: 1

  # Specify if the pipeline should run with a trimer UMI on the tso
  tso_present: 1

  # Specify if a split prefix of index is needed for running minimap2 if the index is large
  minimap2_splitprefix: 0

  # Threshold to remove UMI errors
  error_removal: 1

  # mclumi options

  mclumi:

    editdistance: 9

    memory: 100G


  job_options: -t 48:00:00

**4.** Next we will run the pipleine::

   tallytrin count make full -v5 --no-cluster

This ``--no-cluster`` will run the pipeline locally if you do not have access to a cluster. Alternatively if you have a
cluster remove the ``--no-cluster`` option and the pipleine will distribute your jobs accross the cluster.
