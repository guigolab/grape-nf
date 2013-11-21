Blueprint RNAseq pipeline
=========================

The Blueprint RNAseq pipeline contains all the steps perfomed for the downstream analyses on RNAseq data within the Blueprint project.

The Blueprint RNAseq pipeline is released under GPL. All the programs used in the pipeline are also licensed under GPL except for the GEM binaries 
that are distributed under the `GEM Non commercial binary license <http://algorithms.cnag.cat/wiki/GEM:Non_commercial_binary_license>`_.

Installation
------------

Before installing the Blueprint RNAseq pipeline make sure you have the following software installed in your machine:

- git
- make
- gcc
- python 2.7

You will also need a working internet connection to download and install the pipeline.

Pipeline buildout
-----------------

The first step to install the pipeline is cloning the repository::

    git clone https://github.com/emi80/bp.pipeline.git

This will create a folder in your working directory called ``bp.pipeline``. You should use it as your base folder for running the pipeline. If you want to clone to a different folder you can do::

    git clone https://github.com/emi80/bp.pipeline.git PIPELINE_FOLDER

Once the cloning step completes you can perform the buildout::

    cd PIPELINE_FOLDER
    ./buildout.sh

If your system does not come with python 2.7 installed but you can provide the path to a python 2.7 installation the buildout command can be run like::

    ./buildout.sh /path/to/python2.7

Running the pipeline
--------------------

The software requirements for running the pipeline are:

- bash 4
- awk
- bedtools
- R

An initial setup step is required to produce the following reference files:

#. GEM index of the genome
#. GEM index of the transcriptome
#. fasta index of the genome

To carry out this step you need to activate the python virtualenv from the PIPELINE_FOLDER with ``. bin/activate`` and then run the following commands::

    # create the genome index
    gemtools index -t THREADS -i GENOME

    # create the transcriptome index
    gemtools t-index -t THREADS -i GENOME_INDEX -a ANNOTATION

    # create the fasta index
    samtools faidx GENOME

Once the previous steps are done you will need to prepare the input files. A certain folder structure has to be followed to make the pipeline working. Assuming that you have several samples to be processed you will need the following structure::

    PIPELINE_FOLDER
        <base_files>
        sampleA
            sampleA_1.fastq.gz
            sampleA_2.fastq.gz
        sampleB
            sampleB_1.fastq.gz
            sampleB_2.fastq.gz
        ...
        sampleX
            sampleX_1.fastq.gz
            sampleX_2.fastq.gz

To run the pipeline for a sample you then do::

    cd sampleA
    ../blueprint.pipeline.sh -i sampleA_1.fastq.gz -g GENOME -a ANNOTATION --paired-end [other options]

Run ``../blueprint.pipeline.sh --help`` from the same folder to get the usage message.

Results
-------

Using the previous example and running the pipeline with the default options you will find the follwoing files within you sample folder (e.g. sampleA)::

    PIPELINE_FOLDER/
        <base_files>
        sampleA/
            sampleA_1.fastq.gz
            sampleA_2.fastq.gz
            [sampleA.bedgraph]
            sampleA.bigwig
            sampleA.bigwig.md5
            sampleA_contigs.bed
            sampleA_contigs.bed.md5
            [sampleA.junctions]
            sampleA_m4_n10.bam
            sampleA_m4_n10.bam.bai
            sampleA_m4_n10.bam.bai.md5
            sampleA_m4_n10.bam.md5
            sampleA_m4_n10.map.gz
            sampleA_m4_n10.map.gz.md5
            sampleA_m4_n10.stats
            sampleA_m4_n10.stats.md5
            [sampleA_m4_n10_uniq.bam]
            sampleA.map.gz
            sampleA.map.gz.md5
            stats/
                <stats files>
        quantifications
            sampleA
                sampleA_distinct_exon_with_rpkm.gff
                sampleA_distinct_exon_with_rpkm.gff.md5
                sampleA_flux_profile.log
                sampleA_flux_quantification.log
                sampleA_gene_with_rpkm.gff
                sampleA_gene_with_rpkm.gff.md5
                sampleA.gtf
                sampleA.gtf.md5
                [sampleA_intron.gtf]
                [sampleA_intron.gtf.md5]
                [sampleA_junction.gtf]
                [sampleA_junction.gtf.md5]
                sampleA.par
                sampleA.profile
                sampleA_sort_annotation.log
                sampleA_transcript.gtf
                sampleA_transcript.gtf.md5

The files between brackets could be absent in case a temporary folder has been used in the pipeline run. The ``junction`` and ``intron`` files are created only if the --count-elements parameter contains them. Please refer to the command help for further details.

If the input data is stranded two bigwig files will be present (one for each strand) and will look like::

    [sampleA_m4_n10_1rev.bam]
    [sampleA.plusRaw.bedgraph]
    sampleA.plusRaw.bigwig
    sampleA.plusRaw.bigwig.md5
    [sampleA.minusRaw.bedgraph]
    sampleA.minusRaw.bigwig
    sampleA.minusRaw.bigwig.md5
