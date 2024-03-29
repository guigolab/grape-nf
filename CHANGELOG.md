# GRAPE-nf Changelog

## Version 1.3.0

* Reduce memory usage in RSEM quantification process [78fd890]
* Fix regression - allow usage of precomputed genome index [4a885e3]
* Minor changes to main pipeline script [a41b5dc]
* Fix issue with contig step on unstranded data [87d7db7]
* Fix file issue in fastaIndex process [36f3469]
* Update memory for mapping [e00579a]
* Fix file issues in steps for mapping and quantification indices [574a115]
* Fix syntax errors in rsem module [a7b549e]
* Refactor modules and workflows [df1731f]

## Version 1.2.1 

* Update numeric type for memory in mapping module - fix regression

## Version 1.2 

* Add CFF file for citation - resolve #75
* Update ci script and GH actions workflow
* Skip QC downstream steps for bams with no alignments - fix #76
* Add sample information in process tags - resolve #73
* Correct reported number of sequencing runs for non-unique ids - fix #74
* Update trace file path in main config
* DSL2 refactoring - #78 #62

## Version 1.1.4

* Update cross key to get bamFilesToTranscriptome and bamStatsFiles - fixes #72

## Version 1.1.3

* Handle other == NA in infer experiment script - fix #71

## Version 1.1.2

* Update numeric type for mapping memory - fixes #70

## Version 1.1.1

* Update GitHub workflow
* Set memory limit to STAR index process - fix #66
* Merge pull request #61 from guigolab/dependabot/bundler/docs/addressable-2.8.0
* Bump addressable from 2.4.0 to 2.8.0 in /docs

## Version 1.1.0

* Fix issue with `sambamba` repo checkout in Dockerfile
* Change quantification process view - #59
* Add output files section in readme - #58
* Switch CI to GitHub Actions


## Version 1.0.0

First IHEC release

* Whole pipeline refactoring - #51
* Fix wrong merging information - #55
* Add ability to use regex for reference prefix in `bigwig` step - #17

Additional improvements and fixes:

* Use 2/3 of the process memory in `mapping`
* Increase `quantification` memory in IHEC resources config file
* Update `sambamba` version to `0.7.1`
* Update `bamstats` version to `0.3.4`
* Enable autoMounts for Singularity by default
* Add script to prepare input `FASTQ` files from `BAM` files
* Fix issue in field order when reading input `BAMs`
* Add `gtf` and `gff` to supported file types in CI script
* Sort transcriptome `BAMs` before merge
* Remove `test.sh`, update CI tests and TravisCI configuration file
* Fix workflow for input `BAMs` and add test profile
* Add threshold option to infer experiment script and add corresponding pipeline parameter
* Fix read group attributes in `mapping`

IHEC specific changes:

* Move all IHEC documentation to `ihec-setup.md`, update and clean up
* Clean up pipeline config and use single IHEC image for the IHEC profile
* Add script to pre-fetch the IHEC Singularity image
* Add the IHEC profile to TravisCI config
* Update test data and checksum, add checksums for the IHEC profile

## Version 0.4.1

* Fix wrong comment in config file braking pipeline execution
* Update minimum Nextflow version badge in readme to 0.30.2

## Version 0.4.0

* Add mark-duplicates step - #40.
* Convert existing samtools commands in templates to corresponding available sambamba one - #41.
* Add bamStats process - #42
* Add process and logic to download input files - #43
* Support compressed reference files - #44
* Replace groupBy related logic in pipeline script - #46

Additional fixes and improvements

* Add memory limit to sambamba sort command in RSEM quantification step
* Add sort buffer size option in sambamba markdup
* Update container for fetch process
* Update IHEC profile BAM sorting tool
* Update Nextflow and TravisCI config files
* Update gem mapping templates for new version of samtools
* Add nextflow config profiles for IHEC and update TravisCI configuration
* Clean-up into calls using closure syntax
* Set halfCpus as integer and fix it when task.cpus is 1
* Add Dockerfiles
* Update readme and add description of pipeline parameters
* Add Singularity section to readme with tips and troubleshoot
* Update tool versions in readme
* Bump sambamba version to 0.7.0
* Bump tool versions for bamstats, samtools and RSeQC

## Version 0.3.0

* Fix #21 - sort index file lines when grouping paired `FASTQ` files
* Disable Docker by default and update test script - resolve #22
* Update test profiles and main config - address #25
* Move process resource to external config file - #24
* Modified quantification step by removing the temporary sorted bam file. By doing this we can pipe an uncompressed sam file directly to rsem   - #26 thanks @karl616
* Fix issue with modified quantification templates - introduced by #26
* Use withName process selector in configuration files - close #39
* Update test script, configuration profiles and checksums - #30
* Add option to enforce strandness - first implementation of #11
* Fix bedgraph filename error in contig template - close #27
* Add tool versions to readme - close #18

Additional fixes and imporvements:

* Change `README` markup to Markdown
* Update CI configuration
* Improve `RGCRG` bigwig templates by removing temporary bam files
* Fix missing `NH:1` filter in contigs antisense template
* Move to Travis CI
* Fix issue with test script using wrong pipeline.db when running outside pipeline base dir


## Version 0.2.1

* Decreasing disk stress and processing time of the contig step by replacing temporary bam files by pipes. The drawback is that it is harder to control the maximum number of used cores.
* Add options for docker so that created files are no longer owned by root

## Version 0.2.0

* Added ability to resolve relative paths in the index file against its current location
* Change default profile to `starrsem`
* Set bamSort to `samtools` in default profile
* Fixed Empty lines in the index file result in a `ArrayIndexOutOfBoundsException` #7
* Update documentation - resolve #6
* Add genomeIndex parameter - fix #9
* Add samtools index to mergeBam command template
* Add mapping templates for `STAR 2.5`
* Fix input BAMs considered as paired-end - resolve #14
* Fix error when using both Genome and Transcriptome bams as inputs - fix #13
* Enable trace by default in nextflow config

## Version 0.1.0

First version
