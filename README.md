# Grape

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.23.1-brightgreen.svg)](http://nextflow.io)
[![Build Status](https://travis-ci.org/guigolab/grape-nf.svg?branch=master)](https://travis-ci.org/guigolab/grape-nf)

Grape provides an extensive pipeline for RNA-Seq analyses. It allows the creation of an automated and integrated workflow to manage and analyse RNA-Seq data.

It uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.

## Requirements

- Unix-like operationg system (Linux, MacOS, etc)
- Java 8
- [Docker](https://www.docker.com/) or [Singularity](http://singularity.lbl.gov) engine

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

## Pipeline help

The usage message can be seen using the `--help` pipeline option. 

```
G R A P E ~ RNA-seq Pipeline
----------------------------
Run the GRAPE RNA-seq pipeline on a set of data.

Usage:
    grape-pipeline.nf --index INDEX_FILE --genome GENOME_FILE --annotation ANNOTATION_FILE [OPTION]...

Options:
    --help                              Show this message and exit.
    --index INDEX_FILE                  Index file.
    --genome GENOME_FILE                Reference genome file(s).
    --annotation ANNOTAION_FILE         Reference gene annotation file(s).
    --steps STEP[,STEP]...              The steps to be executed within the pipeline run. Possible values: "mapping", "bigwig", "contig", "quantification". Default: all
    --max-mismatches THRESHOLD          Set maps with more than THRESHOLD error events to unmapped. Default "4".
    --max-multimaps THRESHOLD           Set multi-maps with more than THRESHOLD mappings to unmapped. Default "10".
    --bam-sort METHOD                   Specify the method used for sorting the genome BAM file.
    --add-xs                            Add the XS field required by Cufflinks/Stringtie to the genome BAM file.

SAM read group options:
    --rg-platform PLATFORM              Platform/technology used to produce the reads for the BAM @RG tag.
    --rg-library LIBRARY                Sequencing library name for the BAM @RG tag.
    --rg-center-name CENTER_NAME        Name of sequencing center that produced the reads for the BAM @RG tag.
    --rg-desc DESCRIPTION               Description for the BAM @RG tag.
```

## Configuring the pipeline

### Executors

Nextflow provides different `Executors` to run processes on the local machine, on a computational grid or the cloud without any change to the actual code.
By default a local executor is used, but it can be changed by using [Nextflow executors](http://www.nextflow.io/docs/latest/executor.html).

For example, to run the pipeline in a computational cluster using Sun Grid Engine you can set up a `nextflow.config` file in your current working directory with something like:

```
process {
    executor = 'sge'
    queue    = 'my-queue'
    penv     = 'smp'
}
```

### Input data

The pipeline needs as an input a tab separated file containing containing information about the `FASTQ` files to be processed. The
needed columns in order are:

||||
|-|-|-|
1 | `sample` | the sample identifier, used to merge bam files in case multiple runs for the same sample are present
2 | `id`     | the run identifier (e.g. `test1`)
3 | `path`   | the path to the fastq file
4 | `type`   | the type (e.g. `fastq`)
5 | `view`   | an attribute that specifies the content of the file (e.g. `FqRd` for single-end data or `FqRd1`/`FqRd2` for paired-end data)

**NOTE**: Fastq files from paired-end data will be grouped together by `id`.

Here is an example from the test run:

```
sample1  test1   data/test1_1.fastq.gz   fastq   FqRd1
sample1  test1   data/test1_2.fastq.gz   fastq   FqRd2
```

Sample and id can be the same in case you don't have/know sample identifiers:

```
sample1  test1   data/test1_1.fastq.gz   fastq   FqRd1
sample1  test1   data/test1_2.fastq.gz   fastq   FqRd2
```

### Software

The default Grape configuration uses Docker or Singularity to provision the programs needed for the execution. Pre-built Grape containers are publicly available at the [Grape page in Docker Hub](https://hub.docker.com/r/grape/).

Nextflow also supports [Environment Modules](http://www.nextflow.io/docs/latest/process.html?#module). Creating a working configuration for Grape using modules is not yet straightforward. If you need to use modules please contact us directly.

## Pipeline profiles

The Grape pipeline can be run using different configuration profiles. The profiles essentially allow the user to run the analyses using
different tools and configurations. To specify a profile you can use the [`-profiles` Nextflow option](http://www.nextflow.io/docs/latest/config.html#config-profiles).

The following profiles are available at present:


profile | description
|-|-|
 `gemflux`  | uses `GEMtools` for mapping pipeline and `Flux Capacitor` for isoform expression quantification
 `starrsem` | uses `STAR` for mapping and bigwig and `RSEM` for isoform expression quantification
 `starflux` | uses `STAR` for mapping and `Flux Capacitor` for isoform expression quantification

The default profile uses `STAR` and `RSEM` and set the `--bam-sort` option to `samtools`.

## Run the pipeline

Here is a simple example of how you can run the pipeline:

```
nextflow -bg run grape-nf --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG -resume > pipeline.log
```

By default the pipeline execution will stop as far as one of the processes fails. To change this behaviour you can use the [errorStrategy directive](http://www.nextflow.io/docs/latest/process.html#errorstrategy) of Nextflow processes. You can also specify it on command line. For example, to ignore errors and keep processing you can use `-process.errorStrategy=ignore`.

It is also possible to run a subset of pipeline steps using the option ``--steps``. For example, the following command will only run the ``mapping`` and ``quantification`` steps:

```
nextflow -bg run grape-nf --steps mapping,quantification --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG > pipeline.log
```

## Pipeline results

The pipeline compiles a list of output files into the file `pipeline.db` inside the current working folder. This format of this file is similar to the input with a few more columns:

||||
|-|-|-|
1 | `sample` | the sample identifier, used to merge bam files in case multiple runs for the same sample are present
2 | `id`          | the run identifier (e.g. `test1`)
3 | `path`        | the path to the fastq file
4 | `type`        | the type (e.g. `bam`)
5 | `view`        | an attribute that specifies the content of the file (e.g. `GenomeAlignments`)
6 | `readType`    | the input data type (either `Single-End` or `Paired-End`)
7 | `readStrand`  | the inferred exepriment strandednes if any (it can be `NONE` for unstranded data, `SENSE` or `ANTISENSE` for single-end data, `MATE1_SENSE` or `MATE2_SENSE` for paired-end data.)

Here is an example from the test run:

```
sample1   test1   /path/to/results/sample1.contigs.bed    bed      Contigs                Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.isoforms.gtf   gtf      TranscriptAnnotation   Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.plusRaw.bw     bigWig   PlusRawSignal          Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.genes.gff      gtf      GeneAnnotation         Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/test1_m4_n10.bam       bam      GenomeAlignments       Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.minusRaw.bw    bigWig   MinusRawSignal         Paired-End   MATE2_SENSE
```

## Software versions

The pipeline can be run without a container engine by installing the required software on the local system.

The required tool version for the `standard` profile that have been tested so far are the following:

- [STAR v2.4.0j](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.0j)
- [samtools v1.2](https://github.com/samtools/samtools/releases/tag/1.2)
- [sambamba v0.6.8](https://github.com/biod/sambamba/releases/tag/v0.6.8)
- [RSEM v1.2.21](https://github.com/deweylab/RSEM/releases/tag/v1.2.21)
- [RSeQC v2.6.4](http://rseqc.sourceforge.net/)
- [bedtools v2.19.1](https://github.com/arq5x/bedtools2/releases/tag/v2.19.1)
- [KentUtils v308](https://github.com/ucscGenomeBrowser/kent/releases/tag/v308_base)
- [bamstats v0.3.0](https://github.com/guigolab/bamstats/releases/tag/v0.3.0)