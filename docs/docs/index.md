Title:   Home
Summary: Main page for Grape documentation.
Authors: Emilio Palumbo
Date:    June 15, 2015

# Grape

Grape provides an extensive pipeline for RNA-Seq analyses. It allows the creation of an automated and integrated workflow to manage and analyse RNA-Seq data.

It uses [Nextflow] as the execution backend and can be installed and run in a very simple way. Using Nextflow sharing features the pipeline is automatically downloaded from a git repository and run on the local environment. Please check [Nextflow documentation] for more information.

## Installing Nextflow

Nextflow can be installed with the following commands:

```shell
curl -fsSL get.nextflow.io | bash
mv nextflow <folder in your path>
```

Check the installed version:

```shell
nextflow info
```

```
Version: 0.14.2 build 2998
Modified: 25-06-2015 09:29 UTC (11:29 CEST)
System: Linux 2.6.32-504.8.1.el6.x86_64
Runtime: Groovy 2.3.11 on Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19
Encoding: UTF-8 (UTF-8)
```

!!!note
    Nextflow installation is only required once. You might need to [update](#upgrade-nextflow-version) it, though. 

### Upgrade Nextflow version

The pipeline requires Nextflow version 0.14.0 or nigher. Make sure you have the right version by running `nextflow info`. In case you have an older version you can upgrade as follows:

```shell
nextflow -self-update
```

## Setting up the pipeline

The pipeline will be automatically downloaded the first time if it is not found on the local system. To run the pipeline the following command is used:

```shell
nextflow run guigolab/grape-nf ...
```

The command above will not update the pipeline if it is already present in the system. In order to always run the latest version it is possible to download the pipeline in advance and then run it:

```shell
nextflow pull guigolab/grape-nf
nextflow run grape-nf ...
```

!!! note
    Downloading the pipeline in advance also allows to run the it without specifying the repository `owner` (have a look at [Nextflow documentation](http://www.nextflow.io/docs/latest/sharing.html)).


## Running the pipeline

Using Nextflow, the pipeline can be run on the local machine, on a computational grid or the cloud without any change to the actual code. By default a local executor is used, but it can be changed by using [Nextflow executors](http://www.nextflow.io/docs/latest/executor.html).

For example, to run the pipeline in a computational cluster using Sun Grid Engine you will need to set up a Nextflow configuration file containing something like:

```groovy	
process {
    executor = 'sge'
    queue    = 'myQueue'
    penv     = 'mpi'
}
```

### Creating a working directory

Move to the folder where you want to create the base directory for the pipeline and run the following command:

```shell
mkdir -p <folder name> && cd <folder name>
```

### Reference files

Only genome FASTA file and annotation file are required. Reference files can be specified using absolute paths or relative paths to a local folder where links to the files can be created:

```shell
mkdir refs
ln -s <path to the genome file> refs
ln -s <path to the annotation file> refs
```

### Getting the help

You can get the pipeline help by using the following command:

```shell
nextflow run grape-nf --help
```

```
N E X T F L O W  ~  version 0.14.0
Launching 'guigolab/grape-nf' - revision: 21c5e84cf8 [master]

G R A P E ~ RNA-seq Pipeline
----------------------------
Run the GRAPE RNA-seq pipeline on a set of data.

Usage: 
    grape-pipeline.nf -i INDEX_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]...

Options:
    --help                              Show this message and exit.
    --index INDEX_FILE                  Index file.
    --genome GENOME_FILE                Reference genome file(s).
    --annotation ANNOTAION_FILE         Reference gene annotation file(s).
    --steps STEP[,STEP]...              The steps to be executed within the pipeline run. Possible values: "mapping", "bigwig", "contig", "quantification". Default: all
    --error-strategy ERROR_STRATEGY     Specify how an error condition is managed by the pipeline processes. Possible values: ignore, retry
                                        Default: the entire pipeline  terminates if a process returns an error status.
    --max-mismatches THRESHOLD          Set maps with more than THRESHOLD error events to unmapped. Default "4".
    --max-multimaps THRESHOLD           Set multi-maps with more than THRESHOLD mappings to unmapped. Default "10".

SAM read group options:
    --rg-platform PLATFORM              Platform/technology used to produce the reads for the BAM @RG tag.
    --rg-library LIBRARY                Sequencing library name for the BAM @RG tag.
    --rg-center-name CENTER_NAME        Name of sequencing center that produced the reads for the BAM @RG tag.
    --rg-desc DESCRIPTION               Description for the BAM @RG tag.
```

### Input format

The pipeline needs a tab separated file as an input. This file should contain information about the FASTQ files to be processed. The columns needed in order are:

|
--- | ---
``sample`` | the sample identifier, used to merge bam files in case multiple runs for the same sample are present
``id``     | the run identifier (e.g. labExpId)
``path``   | the path to the fastq file
``type``   | the type (e.g. fastq)
``view``   | an attribute that specifies the content of the file (e.g. FastqRd1)

Here is an example:

```
sample1  test1   data/test1_1.fastq.gz   fastq   FastqRd1
sample1  test1   data/test1_2.fastq.gz   fastq   FastqRd2
```

Sample and id can be the same in case you don't have/know sample identifiers:

```
test  test   data/test.fastq.gz   fastq   FastqRd
```

Also bam files can be specified in the index, with or without fastqs:

```
sample1  test1   data/test1_1.fastq.gz   fastq   FastqRd1
sample1  test1   data/test1_2.fastq.gz   fastq   FastqRd2
sample2  test2   data/test2.bam  bam     Alignment
```

In this case the bam file will skip the mapping process and will be sent to the subsequent processes.

!!! warning
    Bam and fastq files should not refer to the same sample unless you want to merge them!


### Run the pipeline

Here is a simple example of how you can run the pipeline:

```shell
nextflow -bg run grape-nf --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG -resume 2>&1 > pipeline.log
```

By default the pipeline execution will stop as far as one of the processes fails. To change this behaviour you can use the [errorStrategy directive](http://www.nextflow.io/docs/latest/process.html#errorstrategy) of Nextflow processes. You can also specify it on command line. For example, to ignore errors and keep processing you can use ``-process.errorStrategy=ignore``.

It is also possible to run a subset of pipeline steps using the option ``--steps``. For example, the following command will only run the ``mappping`` and ``quantification`` steps:

```shell
nextflow -bg run grape-nf --steps mapping,quantification --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG -resume 2>&1 > pipeline.log
```

### Stop the pipeline

To stop a running pipeline just run the following command from within the pipeline base directory:

```shell    
kill $(cat .nextflow.pid)
```

!!! note
    If you run multiple pipelines within the same folder (e.g. for different genders), please use the [NXF_PID_FILE](http://www.nextflow.io/docs/latest/config.html#environment-variables) environment variable.

### Job monitoring

Nextflow runs all processes in an isolated directory under the pipeline working folder (by default ``./work``). Each process is configured and run by means of several files contained in the process folder. Among those files some can be worth noting:

|
--- | ---
``.command.env`` | the process environment
``.command.out`` | the process standard output
``.command.err`` | the process standard error
``.command.log`` | when run on a compute cluster, the process log output from the job execution
``.command.run`` | the script submitted to the cluster (also contains the header with cluster directives)
``.command.sh``  | the actual command
``.exitcode``    | the exit code of the command

A process can then be easily monitored by inspecting the process folder. Each process is uniquely represented by a hash string nternally computed by Nextflow using comands and inputs. To inspect a process folder just look for Nextflow submission messages in the pipeline log file, which look like the following:

```
...
[b5/0e02e9] Submitted process > index (1)
...
```

The string bewtween square brackets represents the prefix of the relative path to the process folder starting from the ``work`` directory. So to inspect the process workiong folder for the `index (1)` process above:

```shell    
find work/b5 -name '0e02e9*' -exec ls -a {} \+
```

```    
.  ..  .command.begin  .command.env  .command.out  .command.run  .command.sh  .command.val  .exitcode  genome_index.gem  genome_index.log  hg38_AXM.fa
```

## Pipeline profiles

The Grape pipeline can be run using different configuration profiles. The profiles essentially allow the user to run the analyses using different tools. The following profiles are available at present:

|
--- | ---
``gemflux``  | uses the GEMtools for mapping pipeline and the Flux Capacitor for isoform expression quantification
``starrsem`` | uses the STAR for mapping and bigwig and the RSEM program for isoform expression quantification
``starflux`` | uses the STAR for mapping and the Flux Capacitor for isoform expression quantification

The ``gemflux`` profile is used as the default profile if no profile is specified. To specify a profile you can use the `-profiles` Nextflow [option](http://www.nextflow.io/docs/latest/config.html#config-profiles). For example, to use STAR and RSEM just run the following command:

```shell
nextflow -bg run grape-nf -profile starrsem --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG -resume 2>&1 > pipeline.log
```

<!-- Links -->
[Nextflow]: http://nextflow.io
[Nextflow documentation]: http://www.nextflow.io/docs/latest/index.html