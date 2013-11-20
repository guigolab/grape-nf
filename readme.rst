Blueprint RNAseq pipeline
=========================

Installation
------------

Before installing Blueprint RNAseq pipeline make sure you have the following software installed in the machine:

- git
- make
- gcc
- python
- bedtools

You will also need a working internet connection to download and install the pipeline.

Pipeline buildout
-----------------

The first step to install the pipeline is cloning the repository::

    git clone https://github.com/emi80/bp.pipeline.git

This will create a folder in your working directory called ``bp.pipeline``. You should use it as your base folder for running the pipeline. If you want to clone to a different folder you can do::

    git clone https://github.com/emi80/bp.pipeline.git PIPELINE_FOLDER

Once the cloning step completes you can do the buildout::

    cd PIPELINE_FOLDER
    ./buildout.sh

