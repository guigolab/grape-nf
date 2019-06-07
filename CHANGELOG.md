# GRAPE-nf Changelog

## Version 0.3.0

- Fix #21 - sort index file lines when grouping paired FASTQ files
- Disable Docker by default and update test script - resolve #22
- Update test profiles and main config - address #25
- Move process resource to external config file - #24
- Modified quantification step by removing the temporary sorted bam file. By doing this we can pipe an uncompressed sam file directly to rsem   - #26 thanks @karl616
- Fix issue with modified quantification templates - introduced by #26
- Use withName process selector in configuration files - close #39
- Update test script, configuration profiles and checksums - #30
- Add option to enforce strandness - first implementation of #11
- Fix bedgraph filename error in contig template - close #27
- Add tool versions to readme - close #18

Additional fixes and imporvements:

- Change README markup to Markdown
- Update CI configuration
- Improve RGCRG bigwig templates by removing temporary bam files
- Fix missing NH:1 filter in contigs antisense template
- Move to Travis CI
- Fix issue with test script using wrong pipeline.db when running outside pipeline base dir


## Version 0.2.1

- Decreasing disk stress and processing time of the contig step by replacing temporary bam files by pipes. The drawback is that it is harder to control the maximum number of used cores.
- Add options for docker so that created files are no longer owned by root

## Version 0.2.0

- Added ability to resolve relative paths in the index file against its current location
- Change default profile to `starrsem`
- Set bamSort to `samtools` in default profile
- Fixed Empty lines in the index file result in a `ArrayIndexOutOfBoundsException` #7
- Update documentation - resolve #6
- Add genomeIndex parameter - fix #9
- Add samtools index to mergeBam command template
- Add mapping templates for `STAR 2.5`
- Fix input BAMs considered as paired-end - resolve #14
- Fix error when using both Genome and Transcriptome bams as inputs - fix #13
- Enable trace by default in nextflow config

## Version 0.1.0

First version
