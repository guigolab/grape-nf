# IHEC RNA-Seq pipleine

The pipeline developed by [Guigo Lab](https://github.com/guigolab/grape-nf) has been adopted for the IHEC integrative analysis. 
This document describes setting up and testing the pipleine to run on a single machine (or submitted as a batch job to run on a single compute node of a cluster). 

# Setting up and running IHEC RNA-Seq pipeline

## Get nextflow workflow manager
    
    wget https://github.com/nextflow-io/nextflow/releases/download/v19.04.0/nextflow-19.04.0-all -O nextflow
    chmod +x ./nextflow

## Get tarball for Guigo lab RNA-Seq IHEC pipleine

    curl -sL https://github.com/guigolab/grape-nf/archive/IHEC.tar.gz | tar xz
    
## Run initial pipeline tests

Create a pipeline working directory and move to it:

    mkdir pipeline-test && cd pipeline-test

Run the basic testsuite with:

    ../nextflow run ../grape-nf-IHEC -with-singularity 
   
The first time Nextflow will pull and cache all the required Singularity images. By default images will be saved inside the pipeline `work` folder. A different (e.g. shared) location can be specified by either setting the `NXF_SINGULARITY_CACHEDIR` environment variable or creating a file called `nextflow.config` in the current working folder of your pipeline and including the following snippet (replace `<PATH/TO/SINGUALRITY/CACHE/DIR>` with the actual path containing the cached images):

```groovy
singularity {
  cacheDir = "<PATH/TO/SINGUALRITY/CACHE/DIR>"
}
```

Check that the run was successful by looking at status in `trace.txt`

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

The following error

    FATAL:   Unable to pull docker://grapenf/bamstats:bamstats-0.3.2: conveyor failed to get: no descriptor found for reference "cddb8cf488abc7102b1efd0dec0448cd9377ad09606f565013b52251ef9ea1dd"

is a transient connection error that may happen when pulling Singualrity images from the Docker Hub. To solve this, just run the pipeline again until it completes without errors.

## Run the IHEC testsuite defined in the `ihec` nextflow profile:

The pipeline can be run on the MCF10A IHEC RNA-Seq dataset specified in the `ihec` configuration profile. Please note that processing this data requires a considerable amount of [computational resources](https://github.com/guigolab/grape-nf/blob/IHEC/ihec-resources.config). Move back to the initial folder and create a new working directory for the pipeline run:

    cd ..
    mkdir IHEC-pipeline-test && cd IHEC-pipeline-test

Copy the `nextflow.config` file from the initial test folder or create one as needed.

Use the following command to run the IHEC testsuite:

    ../nextflow run ../grape-nf -profile ihec -with-singularity