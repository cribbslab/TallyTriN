.. _getting_started-Installation:


============
Installation
============

The following sections describe how to install the tallytrin workflows. 

.. _getting_started-Conda:

TallyTriN Installation
----------------------

The our preffered method of installation is using conda/mamba. If you dont have conda installed then
please install conda using `miniconda <https://conda.io/miniconda.html>`_ or `anaconda <https://www.anaconda.com/download/#macos>`_.

To install tallytrin::

    conda install mamba -c conda-forge

    # Clone the TallyTriN repo
    git clone https://github.com/cribbslab/TallyTriN.git
    # Use mamba to install environment
    mamba env update --file conda/environment.yml

    conda activate tallytrin
    

    # Install TallyTrinN
    python setup.py install

.. _getting_started-Automated:


Access libdrmaa shared library
------------------------------

You may also need access to the libdrmaa.so.1.0 C library, which can often be installed as part of the
libdrmaa-dev package on most Unixes. Once you have installed that, you may need to tell DRMAA Python
where it is installed by setting the DRMAA_LIBRARY_PATH environment variable, if it is not installed
in a location that Python usually looks for libraries.

In order to set this correctly every time please add the following line to your bashrc, but set the library
path to the location of the libdrmaa.so.1.0::

  export DRMAA_LIBRARY_PATH=/usr/lib/libdrmaa.so.1.0



.. _conda: https://conda.io
