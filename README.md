# TallyNNN



Overview
========
Long-read sequencing has become an increasingly popular tool for RNA sequencing that provides unprecedented insight into isoform, translocation and variant calling analysis.

Tally-NNN is a collection of bulk and single-cell workflows that utlise Unique Molecular Identifiers (UMIs) synthesised using trimer blocks of nucleosides.

Workflows
=========
Included within this repo are three workflows:

* pipeline_count - a bulk RNA-seq workflow that facilitates the analysis of UMIs synthesised using trimer nucleoside blocks added to the 5' and 3' of the RNA molecule.
* pipeline_fusion - a bulk analysis pipeline for the analysis of fusion transcripts. It implements a strategy for removing chimeric artefacts and quantifies real genomic translocation at the RNA-seq level.
* pipeline_singlecell - a workflow that facilitates the analysis of scCOLOR-seq nanopore single-cell sequencing data, but with trimers added at the 5' end of the RNA. 

Installation
============

We reccomend installing [miniconda](https://docs.conda.io/en/latest/miniconda.html), then creating
a new environment and install mamba

  ```
  conda install mamba -c conda-forge
  ```
  
Next install the required software using the conda yml file 

  ```
  mamba env update --file conda/environment.yml
  ```

Activate the condda environment

  ```
  conda activate trimer
  ```

Then, you will need to manually install TallyNN and the fork of umi tools. The fork is added as a submodule to this
repo to help you easily install.

  ```
  # Clone the TallyNNN repo
  git clone 
  # Install TallyNNN code
  python setup.py install
  ```

Documentation
=============

Further information how you can run the pipelines can be found at [read the docs](https://tallynnn.readthedocs.io/en/latest/)



Usage
=====

Run the ``tallynn --help`` command to see what workflows are available and ``tallynnn count  -help`` to see how to use them.


For example, to generate a configuration file run

   ```
   tallynnn count config
   ```

To set up the configuration file please refer to [read the docs](https://tallynnn.readthedocs.io/en/latest/getting_started/Tutorial.html#modify-the-config-file).

To run the pipeline with all tasks then run
   
   ```
   tallynnn count make full -v5 
   ```

Manuscript
==========

The bioRxiv manuscript accompanying this code can be found here: 


```