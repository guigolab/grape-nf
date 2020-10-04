# IHEC RNA-Seq pipeline

The pipeline developed by [Guigo Lab](https://github.com/guigolab/grape-nf) has been adopted for the IHEC integrative analysis.
This document describes setting up and testing the pipeline to run on a single machine (or submitted as a batch job to run on a single compute node of a HPC cluster environment). To run the IHEC pipeline at scale on a HPC cluster environment taking advantage of a job scheduler please check [Nextflow documentation](https://www.nextflow.io/docs/latest/basic.html#execution-abstraction) or contact the pipeline developer.

## Setting up and running IHEC RNA-Seq pipeline

### Requirements

- Java Runtime 8 or newer
- [Singularity](https://sylabs.io/singularity/) 3.0 or newer

### Get the Nextflow workflow manager

    wget https://github.com/nextflow-io/nextflow/releases/download/v19.04.0/nextflow-19.04.0-all -O nextflow
    chmod +x ./nextflow
    # optional: move the nextflow executable in a folder within your PATH environment variable
    mv ./nextflow <FOLDER/WITHIN/PATH>

### Get tarball for Guigo lab RNA-Seq IHEC pipeline

    mkdir grape-nf-IHEC
    curl -sL https://github.com/IHEC/grape-nf/archive/master.tar.gz | tar xz --strip-components=1 -C grape-nf-IHEC

### Run initial pipeline tests

Create a pipeline working directory and move to it:

    mkdir pipeline-test && cd pipeline-test

A specific pipeline profile with the standard workflow recommended by the IHEC consortium for RNA-seq data is available. The profile is on continuous integration testing and can be used by specifying the option `-profile ihec` on the command line. To run the IHEC profile with a small test dataset type:

    ../nextflow run ../grape-nf-IHEC -profile ihec -with-singularity

The first time Nextflow will pull and cache all the required Singularity image. By default it will be saved inside the pipeline `work` folder. A different (e.g. shared) location can be specified by either setting the `NXF_SINGULARITY_CACHEDIR` environment variable or creating a file called `nextflow.config` in the current working folder of your pipeline and including the following snippet (replace `<PATH/TO/SINGULARITY/CACHE/DIR>` with the actual path containing the cached images):

    singularity {
        cacheDir = "<PATH/TO/SINGULARITY/CACHE/DIR>"
    }

**NOTE**: As an option, the script `../grape-nf-IHEC/scripts/singularity-prefetch.sh` can be used to pre-fetch the singularity image in advance. Run `../grape-nf-IHEC/scripts/singularity-prefetch.sh -h` for usage.

Check that the run was successful by looking at the `status` field in `trace.txt`:

    $ cat trace.txt
    task_id hash    native_id       name    status  exit    submit  duration        realtime        %cpu    peak_rss        peak_vmem       rchar   wchar
    1       05/6d4ea8       183145  fastaIndex (genome-SAMtools-1.3.1)      COMPLETED       0       2019-06-22 12:10:31.913 532ms   58ms    99.1%   0       0       1.3 MB  249 B
    3       16/de6889       183195  txIndex (genome-RSEM-1.2.21)    COMPLETED       0       2019-06-22 12:10:31.961 706ms   301ms   85.8%   9 MB    23 MB   3.5 MB  4.9 MB
    2       b1/928996       183164  index (genome-STAR-2.4.0j)      COMPLETED       0       2019-06-22 12:10:31.941 12s     11.6s   70.0%   2.7 GB  2.7 GB  26.9 MB 1.5 GB
    4       4e/639dbd       190888  mapping (test1-STAR-2.4.0j)     COMPLETED       0       2019-06-22 12:10:44.090 3s      2.6s    101.0%  1.1 GB  1.7 GB  1.5 GB  551.2 KB
    5       6e/c564db       193327  inferExp (test1-RSeQC-2.6.4)    COMPLETED       0       2019-06-22 12:10:47.230 1.6s    1.1s    98.4%   3.6 MB  14.5 MB 5.2 MB  101.7 KB
    6       b3/8f9b9f       1551    bamStats (test1-bamstats-0.3.2) COMPLETED       0       2019-06-22 12:10:53.127 557ms   234ms   131.1%  14.3 MB 111.9 MB        645.8 KB        7.1 KB
    8       6a/2b83cf       2258    bigwig (test1-STAR-2.4.0j)      COMPLETED       0       2019-06-22 12:10:54.021 705ms   259ms   77.1%   44.5 MB 69.5 MB 929.8 KB        394.5 KB
    7       ef/bc41f1       4382    contig (test1-RGCRG-0.1)        COMPLETED       0       2019-06-22 12:10:57.139 1.5s    1s      108.2%  1 MB    6 MB    6.3 MB  2.5 MB
    9       7d/dbbd06       4888    quantification (test1-RSEM-1.2.21)      COMPLETED       0       2019-06-22 12:10:57.571 8.8s    8.4s    137.7%  73 MB   1 GB    75 MB   63.8 MB

### Troubleshooting

The following error:

    FATAL:   Unable to pull docker://grapenf/<IMAGE>: conveyor failed to get: no descriptor found for reference "<HASH>"

is a transient connection error that may happen when pulling Singularity images from the Docker Hub. To solve this, just run the pipeline again until it completes without errors or use the script described in the [section above](#run-initial-pipeline-tests) to pre-fetch the image.

## Run the IHEC test dataset

A reference RNA-seq dataset (MCF10A) is made available by the IHEC consortium for testing. Please note that processing this data requires a considerable amount of [computational resources](https://github.com/IHEC/grape-nf/blob/master/config/resources/ihec.config).

**NOTE**: When running the pipeline with real data please make sure your computation environment has enough space for temporary data. Some steps of the pipeline make use of `/tmp` or `$TMPDIR` to store temporary files and need a certain amount of space. In order to specify a custom temporary directory you can set the `TMPDIR` environment variable in your local `nextflow.config` file (please see [Nextflow docs](https://www.nextflow.io/docs/edge/config.html#scope-env) for more information), e.g.:
  ```
  env {
      TMPDIR = "/scratch/tmp"
  }
  ```

Move back to the initial folder and create a new working directory for the pipeline run:

    cd ..
    mkdir IHEC-pipeline-test && cd IHEC-pipeline-test

Copy the `nextflow.config` file from the initial test folder or create one as needed.

Use the following command to run the IHEC test suite:

    ../nextflow run ../grape-nf-IHEC -profile ihec,ihec-data -with-singularity

**NOTE**: the `ihec-data` profile only contains configuration for the required input references and data. Thus it must be used in conjuction with a workflow profile such as `ihec`, as it can be seen in the example command above. Please check [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) for more information about configuration profiles.

Once the data has been processed the following checksums can be used to verify that the pipeline ran fine and produced the correct output files:

    ea4075877567c6b36061370d50fde9e2 *A24901.contigs.bed
    4804ef82c86baf5af551c0ae03204448 *A24901.genes.results
    a111578132a844a34760a9f1ab308a05 *A24901.isoforms.results
    016ef39b20d8d98ce9bf0d789e7225f8 *A24901_stats.json
    3b1489a9452ebfaab2aae9b6e3c47523 *A24901.UniqueMultiple.raw.bw
    2deefb50d097d7f4d66458b40b7e80d3 *A24901.Unique.raw.bw
    5bace482e2c8907e1287d747d4eacb59 *A24901_m4_n10_toGenome.markdup.bam
    c5f698ec8956be3a5fbe3a7325e5c4df *A24901_m4_n10_toTranscriptome.bam

## Run the IHEC pipeline on custom data

In order to run the IHEC pipeline on custom data you need to prepare the input index file following the [Pipeline input](https://github.com/guigolab/grape-nf#pipeline-input) section of the readme.

Create a new folder for the run and move there. If needed, copy the `nextflow.config` file from the initial test folder or create a new one following the steps above. Once the index file is ready, and assuming that you called it `input-files.tsv`, you can use the following command to run the workflow:

    ../nexflow -bg run ../grape-nf-IHEC -profile ihec,ihec-data -with-singularity --index input-files.tsv > pipeline.log
