# GRAPE-nf Changelog

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
