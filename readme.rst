Blueprint RNAseq pipeline
=========================

Installation
------------

Before installing Blueprint RNAseq pipeline make sure you have the following software installed in the machine:

- git
- make
- gcc
- python 2.7
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

If your system does not come with python 2.7 and you ar able to provide the pathe of a python 2.7 binary you can run the buildout in the following way::

    ./buildout.sh /path/to/python2.7

Running the pipeline
--------------------

In order to run the pipeline a setup step is required to produce the following files:

#. GEM index of the genome
#. GEM index of the transcriptome
#. fasta index of the genome

You need to run the following commands::

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
    ../blueprint.pipeline.sh -i sampleA_1.fastq.gz -g GENOME -a ANNOTATION ...other options...

Run ``../blueprint.pipeline.sh --help`` from the same folder to get the usage message.

Results
-------

The results are organized in the following structure. Using the previous example and running the pipeline with the default options you will have something like::

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
                sampleA_intron.gtf
                sampleA_intron.gtf.md5
                sampleA_junction.gtf
                sampleA_junction.gtf.md5
                sampleA.par
                sampleA.profile
                sampleA_sort_annotation.log
                sampleA_transcript.gtf
                sampleA_transcript.gtf.md5

The files in brackets could be absent in case you run the pipeline specifying a temporary folder.

If the input data is stranded two bigwig files will be present and will look like::

    [sampleA_m4_n10_1rev.bam]
    [sampleA.plusRaw.bedgraph]
    sampleA.plusRaw.bigwig
    sampleA.plusRaw.bigwig.md5
    [sampleA.minusRaw.bedgraph]
    sampleA.minusRaw.bigwig
    sampleA.minusRaw.bigwig.md5








