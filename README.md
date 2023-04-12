# TallyTriN



Overview
========
Long-read sequencing has become an increasingly popular tool for RNA sequencing that provides unprecedented insight into isoform, translocation and variant calling analysis.

TallyTriN is a collection of bulk and single-cell workflows that utilize Unique Molecular Identifiers (UMIs) synthesized using trimer blocks of nucleotides.

Workflows
=========
Included within this repo are the follwing workflows:

* pipeline_count - a bulk RNA-seq workflow that facilitates the analysis of UMIs synthesised using trimer nucleoside blocks added to the 5' and 3' of the RNA molecule.
* pipeline_fusion - a bulk analysis pipeline for the analysis of fusion transcripts. It implements a strategy for removing chimeric artefacts and quantifies real genomic translocation at the RNA-seq level.
* pipeline_singlecell - a workflow that facilitates the analysis of scCOLOR-seq nanopore single-cell sequencing data, but with trimers added at the 5' end of the RNA. 
* pipeline_10X - a workflow that facilitates the analysis of 10X sequencing data.

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
  conda activate tallytrin
  ```

Then, you will need to manually install TallyTriN and the fork of umi tools. The fork is added as a submodule to this
repo to help you easily install.


```
# Clone the TallyTriN repo
git clone git@github.com:cribbslab/TallyTriN.git

# Install TallyTriN code
python setup.py install
```

Documentation
=============

Further information how you can run the pipelines can be found at [read the docs](https://tallynnn.readthedocs.io/en/latest/)



Usage
=====

Run the tallytrin --help command to see what workflows are available and tallytrin count -help to see how to use them.

For example, to generate a configuration file run

   ```
   tallytrin count config
   ```

To set up the configuration file, please refer to our [documentation]().

To run the pipeline with all tasks, run

   
   ```
  tallytrin count make full -v5 

   ```

Manuscript
==========

The bioRxiv manuscript accompanying this code can be found here: 

https://www.biorxiv.org/content/10.1101/2023.04.06.535911v1

```
