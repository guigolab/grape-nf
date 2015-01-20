===============
Getting started
===============

Creating a working directory
----------------------------

Move to the folder where you want to create the base directory for the pipeline and run the following command:

.. code-block:: bash

    $ git clone -b nextflow /users/rg/projects/git/bp.pipeline.git <pipeline name>
    $ cd <pipeline name>


Reference files
---------------

Only genome FASTA file and annotation file are required. Reference files can be specified using absolute paths or relative paths to a local folder where links to the files can be created:

.. code-block:: bash

    $ mkdir refs
    $ ln -s <path to the genome file> refs
    $ ln -s <path to the annoation file> refs

.. note:: The genome and the annotation needs to be sorted by **genomic position**.


Installing Nextflow
-------------------
This step will install nextflow for your user. It is only needed once. These are the commands:

.. code-block:: bash

    $ curl -fsSL get.nextflow.io | bash
    $ mv nextflow <folder in your path>
    $ nextflow info
    Version: 0.12.0 build 2580
    Modified: 05-01-2015 15:07 UTC (16:07 CEST)
    System: Linux 2.6.32-504.1.3.el6.x86_64
    Runtime: Groovy 2.3.9 on Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19
    Encoding: UTF-8 (UTF-8)


Pipeline usage
--------------
To get the usage string and the list of options use this command:

.. code-block:: bash

    $ nextflow run blueprint.pipeline.nf --help
    N E X T F L O W  ~  version 0.12.0

    B L U E P R I N T ~ RNA Pipeline
    --------------------------------
    Run the RNAseq pipeline on a set of data.
    
    Usage: 
        blueprint.pipeline.nf -i INDEX_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]...
    
    Options:
        --help                              Show this message and exit.
        --index INDEX_FILE                  Index file.
        --genome GENOME_FILE                Reference genome file(s).
        --annotation ANNOTAION_FILE         Reference gene annotation file(s).
        --chunk-size                        The number of records to be put in each chunk when splitting the input. Default: no split
        --error-strategy                    Specify how an error condition is managed by the pipeline processes. Default: the entire pipeline
                                            terminates if a process returns an error status.
        --max-read-length READ_LENGTH       The maximum read length (used to compute the transcriptomes). Default: "auto".
        --max-mismatches THRESHOLD          Set maps with more than <threshold> error events to unmapped. Default "4".
        --max-multimaps THRESHOLD           Set multi-maps with more than <threshold> mappings to unmapped. Default "10".
        --filter-intron-length THRESHOLD    Filter multimaps preferring ones with intron length > THRESHOLD
        --filter-block-length THRESHOLD     Filter multimaps preferring ones with block length > THRESHOLD
        --filter-level LEVEL                Reduce multimaps using the specified uniqueness level.
        --rg-platform PLATFORM              Platform/technology used to produce the reads for the BAM @RG tag.
        --rg-library LIBRARY                Sequencing library name for the BAM @RG tag.
        --rg-center-name CENTER_NAME        Name of sequencing center that produced the reads for the BAM @RG tag.
        --rg-desc DESCRIPTION               Description for the BAM @RG tag.
        --flux-profile                      Specify whether the Flux Capacitor profile file shoudl be written. Default: "false".
        --count-elements ELEMENTS           A comma separated list of elements to be counted by the Flux Capacitor.
                                            Possible values: INTRONS,SPLICE_JUNCTIONS. Defalut: "none".


Input format
------------

The pipeline needs a tab separated file as one of the inputs. This files should contain information about the FASTQ files to be processed. The columns needed in order are:

==========  ====================================================================================================
``sample``  the sample identifier, used to merge bam files in case multiple runs for the same sample are present
``id``      the run identifier (e.g. labExpId)
``path``    the path to the fastq file
``type``    the type (e.g. fastq)
``view``    an attribute that specifies the content of the file (e.g. FastqRd1)
==========  ====================================================================================================

Here is an example:: 

   sample1  test1   data/test1_1.fastq.gz   fastq   FastqRd1
   sample1  test1   data/test1_2.fastq.gz   fastq   FastqRd2

Sample and id can be the same in case you don't have/know sample identifiers::

   test  test   data/test.fastq.gz   fastq   FastqRd

Run the pipeline
----------------

Here is a simple example of the command to run the pipeline:

.. code-block:: bash

    $ nextflow -bg run blueprint.pipeline.nf --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG -resume 2>&1 > pipeline.log

Stop the pipeline
-----------------

To stop a running pipeline just run the following command:

.. code-block:: bash

    $ kill $(cat .nextflow.pid)

.. note:: If you run multiple pipelines within the same folder (e.g. for different genders), the file ``.nextflow.pid`` will contain only the pid of the last Nextflow pipeline (*while all pipelines will still be running*). In that case you will need to manually look for the pid and kill the corresponding process. Since this approach can be quite confusing it is strongly discouraged for the moment. The problem has been reported and it is likely to be fixed in a next release.

Job monitoring
---------------
Nextflow runs all processes in an isolated directory under the pipeline working folder (by default ``./work``). Each process is configured and run by means of several files contained in the process folder. Among those files some can be worth noting:

================  ======================================================================================
``.command.env``  the process environment
``.command.out``  the process output (merge ``stderr`` into ``stdout`` if not redirected)
``.command.run``  the script submitted to the cluster (also contains the header with cluster directives)
``.command.sh``   the actual command
``.exitcode``     the exit code of the command
================  ======================================================================================

A processe can then be easily monitored by inspecting the process folder. Each process is uniquely represented by a hash string nternally computed by Nextflow using comands and inputs. To inspect a process folder just look for Nextflow submission messages in the pipeline log file, which look like the following::

    ...
    [b5/0e02e9] Submitted process > index (1)
    ...

The string bewtween square brackets represents the prefix of the relative path to the process folder starting from the ``work`` directory. So to inspect the process workiong folder for the `index (1)` process above:

.. code-block:: bash

    $ find work/b5 -name '0e02e9*' -exec ls -a {} \+
    .  ..  .command.begin  .command.env  .command.out  .command.run  .command.sh  .command.val  .exitcode  genome_index.gem  genome_index.log  hg38_AXM.fa


