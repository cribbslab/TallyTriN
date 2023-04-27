.. _manual-main:

========================
TallyTriN documentation!
========================

Long-read sequencing has become an increasingly popular tool for RNA sequencing that provides unprecedented insight into isoform, translocation and variant calling analysis.

TallyTriN is a collection of bulk and single-cell workflows that utilize Unique Molecular Identifiers (UMIs) synthesized using trimer blocks of nucleotides.

Workflows
=========
Included within this repo are the follwing workflows:

* pipeline_count - a bulk RNA-seq workflow that facilitates the analysis of UMIs synthesised using trimer nucleoside blocks added to the 5' and 3' of the RNA molecule.
* pipeline_fusion - a bulk analysis pipeline for the analysis of fusion transcripts. It implements a strategy for removing chimeric artefacts and quantifies real genomic translocation at the RNA-seq level.
* pipeline_singlecell - a workflow that facilitates the analysis of scCOLOR-seq nanopore single-cell sequencing data, but with trimers added at the 5' end of the RNA. 
* pipeline_10X - a workflow that facilitates the analysis of 10X sequencing data.

.. _manual-quick_example:

--------
Citation
--------

The bioRxiv manuscript accompanying this code can be found here:: 
   
   https://www.biorxiv.org/content/10.1101/2023.04.06.535911v1

.. _manual-support:


.. toctree::
   :caption: Getting started
   :name: getting_started
   :maxdepth: 1
   :hidden:

   getting_started/Installation.rst
   getting_started/Tutorial.rst

.. toctree::
   :caption: Project Info
   :name: project-info
   :maxdepth: 1
   :hidden:

   project_info/Licence.rst
