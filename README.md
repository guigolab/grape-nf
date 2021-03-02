# Grape

[![Nextflow version: >=0.30.2](https://img.shields.io/badge/nextflow-%E2%89%A50.30.2-brightgreen.svg)](http://nextflow.io)
[![Singularity version: >=3.x](https://img.shields.io/badge/singularity-3.x-blue.svg)](http://sylabs.io/singularity)
[![CI status](https://github.com/guigolab/grape-nf/workflows/Pipeline%20CI/badge.svg)](https://github.com/guigolab/grape-nf/actions)

Grape provides an extensive pipeline for RNA-Seq analyses. It allows the creation of an automated and integrated workflow to manage and analyse RNA-Seq data.

It uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.

Grape has been adopted for RNA-seq integrative analysis within the IHEC consortium. Check the [IHEC setup document](ihec-setup.md) to run the pipeline following IHEC recommendations.

## Requirements

- Unix-like operationg system (Linux, MacOS, etc)
- Java 8
- [Docker](https://www.docker.com/)  or [Singularity](http://singularity.lbl.gov) engine

## Quickstart

1. Install Nextflow by using the following command:

    ```
    curl -s https://get.nextflow.io | bash
    ```

2. Make a test run:

    ```
    nextflow run guigolab/grape-nf -with-docker
    ```

**NOTE**: the very first time you execute it, it will take a few minutes to download the pipeline from this GitHub repository and the associated Docker images needed to execute the pipeline.

## Pipeline software

The preferred way to run the pipeline is to use Docker or Singularity to provision the programs needed for the execution. Just use the `-with-docker` or `-with-singularity` option in the pipeline command. Pre-built Grape containers are publicly available at the [Grape page in Docker Hub](https://hub.docker.com/r/grape/).

Alternatively, a Conda [environment file](conda/grape-env.yml) is available to create a Conda environment with all the required software. It can either be used to create the environemnt in advance with:

```
conda env create -f conda/grape-env.yml
```

or to delegate Nextflow to prepare the environment using the conda configuration [directive](https://www.nextflow.io/docs/latest/conda.html?highlight=conda#use-conda-environment-files).

### Using Singularity

Singularity is the preferred container engine for running the pipeline in an HPC environment. In order to minimize the amount of issues that could arise we recommend the use of Singularity version 3.0 or higher. At the time of writing of this document the latest Singularity release is 3.2.

#### Image cache dir

The first time you run the pipeline with Singularity it will download the required images from the [Docker Hub](https://hub.docker.com/) and save them in a folder inside the pipeline work dir. You can specify a different location (e.g. a centralized cache) by using the `NXF_SINGULARITY_CACHEDIR` environment variable or by including the following snippet in a file called `nextflow.config` and placing it in the current working folder of your pipeline:

```
singularity {
  cacheDir = "/data/singularity"
}
```

Please check the [Singularity section](https://www.nextflow.io/docs/latest/singularity.html) in Nextflow documentation for more information.

#### Bind mounts

Nextflow expects that data paths are defined system wide, and your Singularity images need to be able to access these paths. Singularity allows paths that do not currently exist within the container to be created and mounted dynamically by specifying them on the command line. For this to work the [user bind control](https://singularity-admindoc.readthedocs.io/en/latest/the_singularity_config_file.html#user-bind-control-boolean-default-yes) option must be set to `true` in the Singularity config file. Nextflow support for this feature is enabled by default for the pipeline, by defining the `singularity.autoMounts = true` setting in the main configuration file.

Starting in version 3.0, Singularity can bind paths to non-existent mount points within the container even in the absence of the “overlay fs” feature, thus supporting architectures running legacy kernel versions (including RHEL6 vintage systems). For older versions of Singularity a kernel supporting the OverlayFS union mount filesystem is required for this functionality to be supported.

Please see [here](https://www.sylabs.io/guides/3.2/user-guide/bind_paths_and_mounts.html) for further instructions on Singularity mounts.

## Pipeline parameters

A usage message is provided and can be seen using the `--help` pipeline option in the command as follows:

```
nextflow run guigolab/grape-nf --help
```

#### `--index INDEX_FILE`

- specifies the path of the file containing the list of input files and the corresponding metadata (see the [next section](#pipeline-input) for more details).

#### `--genome GENOME_FILE`

- sets the location of the input genome `FASTA` file

#### `--annotation ANNOTATION_FILE`

- sets the location of the input `GTF`/`GFF` annotation file

#### `--steps STEP[,STEP]..`

- defines the pipeline steps to be performed

#### `--paired-end`

- specifies that the data is paired-end (to be used whith `BAM` input files)

### Mapping options

#### `--max-mismatches THRESHOLD`

- set a maximum threashold for the number of allowed mismatches

#### `--max-multimaps THRESHOLD`

- set a maximum threashold for the number of allowed multiple mapped reads

#### `--bam-sort METHOD`

- set the sort method of the out `BAM` file

#### `--add-xs`

- add the `SAM` tag `XS` to the output `BAM` file (useful for using the file with tools like Cufflinks or StringTie that use the tag to know the directionality of the split maps)

### Read group options

These options are used to customize the `@RG` header tag of the `BAM` files produced by the mapping step, according to the [`SAM` specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).

#### `--rg-platform PLATFORM`

- set the `PL` attribute

#### `--rg-library LIBRARY`

- set the `LB` attribute

#### `--rg-center-name CENTER_NAME`

- set the `CN` attribute

#### `--rg-desc DESCRIPTION`

- set the `DS` attribute

## Pipeline input

The pipeline reads the paths of the `FASTQ`/`BAM` files to be processed and the corresponding metadata from a `TSV` file (see the `--index` parameter). The file must contains the following columns in order:

||||
|-|-|-|
1 | `sampleID` | the sample identifier, used to merge bam files in case multiple sequencing runs of the same sample are present
2 | `runID`     | the run identifier (e.g. `test1`)
3 | `path`   | the path to the fastq file (it can be absolute or relative to the `TSV` file)
4 | `type`   | the type (e.g. `fastq`)
5 | `view`   | an attribute that specifies the content of the file (e.g. `FqRd` for single-end data or `FqRd1`/`FqRd2` for paired-end data)

**NOTE**: Fastq files from paired-end data will be grouped together by `runID`.

Here is an example from the test run:

```
sample1  test1   data/test1_1.fastq.gz   fastq   FqRd1
sample1  test1   data/test1_2.fastq.gz   fastq   FqRd2
```

Sample and id can be the same in case you don't have/know sample identifiers:

```
run1  run1   data/test1_1.fastq.gz   fastq   FqRd1
run1  run1   data/test1_2.fastq.gz   fastq   FqRd2
```

## Pipeline results

The paths of the resulting output files and the corresponding metadata are stored into the `pipeline.db` file (`TSV` formatted) which sits inside the current working folder. The format of this file is the same as the index file with few more columns:

||||
|-|-|-|
1 | `sampleID` | the sample identifier, used to merge bam files in case multiple runs for the same sample are present
2 | `runID`          | the run identifier (e.g. `test1`)
3 | `path`        | the path to the fastq file
4 | `type`        | the type (e.g. `bam`)
5 | `view`        | an attribute that specifies the content of the file (e.g. `GenomeAlignments`)
6 | `readType`    | the input data type (either `Single-End` or `Paired-End`)
7 | `readStrand`  | the inferred experiment strandedness if any (it can be `NONE` for unstranded data, `SENSE` or `ANTISENSE` for single-end data, `MATE1_SENSE` or `MATE2_SENSE` for paired-end data.)

Here is an example from the test run:

```
sample1   test1   /path/to/results/sample1.contigs.bed    bed      Contigs                     Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.isoforms.gtf   gtf      TranscriptQuantifications   Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.plusRaw.bw     bigWig   PlusRawSignal               Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.genes.gff      gtf      GeneQuantifications         Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/test1_m4_n10.bam       bam      GenomeAlignments            Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.minusRaw.bw    bigWig   MinusRawSignal              Paired-End   MATE2_SENSE
```

### Output files

The pipeline produces several output files during the workflow execution. Many files are to be considered temporary and can be removed once the pipeline completes. The following files are the ones reported in the `pipeline.db` file and are to be considered as the pipeline final output.

#### Alignments to the reference genome

|views|
|-|
|`GenomeAlignments`|

This BAM file contains information on the alignments to the reference genome. It includes all the reads from the FASTQ input. Reads that do not align to the reference are set as unmapped in the bam file. The file can be the product of several steps of the pipeline depending on the given input parameters. It is initially produced by the `mapping` step, then it can be the result of merging of different runs from the same experiment and finally it can run through a marking duplicates process that can eventually remove reads that are marked as duplicates.

#### Alignments to the reference transcriptome

|views|
|-|
|`TranscriptomeAlignments`|

This BAM file contains information on the alignments to the reference transcriptome. It is generally used only for expression abundance estimation, as input in the `quantification` process. The file is generally produced in the `mapping` process and can be the result of merging of different runs from the same experiment.

#### Alignments statistics

|views|
|-|
|`BamStats`|

This JSON file contains alignment statistics computed with the [bamstats](https://github.com/guigolab/bamstats) program. It also reports RNA-Seq quality check metrics agreed within the IHEC consortium.

#### Signal tracks

|views|
|-|
|`RawSignal`|
|`MultipleRawSignal`|
|`MinusRawSignal`|
|`PlusRawSignal`|
|`MultipleMinusRawSignal`|
|`MultiplePlusRawSignal`|

These BigWig files (one or two, depending on the strandedness of the input data) represent the RNA-Seq signal.

#### Contigs

|views|
|-|
|`Contigs`|

This BED file reports RNA-seq contigs computed from the pooled signal tracks.

#### Quantifications

|views|
|-|
|`GeneQuantifications`|
|`TranscriptQuantifications`|

These two files report abundances for genes and transcripts in the processed RNA-seq samples. The format can be either GFF or TSV depending on the tool used to perform the quantification.

## Pipeline configuration

### Executors

Nextflow provides different `Executors` run the processes on the local machine, on a computational cluster or different cloud providers without the need to change the pipeline code.

By default the local executor is used, but it can be changed by using  the [executor](https://www.nextflow.io/docs/latest/config.html#scope-executor) configuration scope.

For example, to run the pipeline in a computational cluster using Sun Grid Engine you can create a `nextflow.config` file in your current working directory with something like:

```
process {
    executor = 'sge'
    queue    = 'my-queue'
    penv     = 'smp'
}
```

### Pipeline profiles

The Grape pipeline can be run using different configuration profiles. The profiles essentially allow the user to run the analyses using
different tools and configurations. To specify a profile you can use the [`-profiles` Nextflow option](http://www.nextflow.io/docs/latest/config.html#config-profiles).

The following profiles are available at present:


profile | description
|-|-|
 `gemflux`  | uses `GEMtools` for mapping pipeline and `Flux Capacitor` for isoform expression quantification
 `starrsem` | uses `STAR` for mapping and bigwig and `RSEM` for isoform expression quantification
 `starflux` | uses `STAR` for mapping and `Flux Capacitor` for isoform expression quantification

The default profile is `starrsem`.

## Run the pipeline

Here is a simple example of how you can run the pipeline:

```
nextflow -bg run grape-nf --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG -resume > pipeline.log
```

By default the pipeline execution will stop as far as one of the processes fails. This behaviour can be changed using the [errorStrategy](http://www.nextflow.io/docs/latest/process.html#errorstrategy) process directive, which can also be specified on the command line. For example, to ignore errors and keep processing you can use:

`-process.errorStrategy=ignore`.

It is also possible to run a subset of the pipeline steps using the option ``--steps``. For example, the following command will only run ``mapping`` and ``quantification``:

```
nextflow -bg run grape-nf --steps mapping,quantification --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG > pipeline.log
```

##  Tools versions

The pipeline can be also run natively by installing the required software on the local system or by using [Environment Modules](http://www.nextflow.io/docs/latest/process.html?#module).

The versions of the tools that have been tested so far with the `standard` pipeline profile are the following:

- [bamstats v0.3.4](https://github.com/guigolab/bamstats/releases/tag/v0.3.4)
- [bedtools v2.19.1](https://github.com/arq5x/bedtools2/releases/tag/v2.19.1)
- [KentUtils v308](https://github.com/ucscGenomeBrowser/kent/releases/tag/v308_base)
- [RSEM v1.2.21](https://github.com/deweylab/RSEM/releases/tag/v1.2.21)
- [RSeQC v2.6.4](http://rseqc.sourceforge.net/)
- [sambamba v0.7.1](https://github.com/biod/sambamba/releases/tag/v0.7.1)
- [samtools v1.3](https://github.com/samtools/samtools/releases/tag/1.3.1)
- [STAR v2.4.0j](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.0j)
