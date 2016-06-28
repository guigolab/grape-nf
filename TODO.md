# Things to do

- [ ] add FastQC
- [ ] add preprocessing steps for input (fastq) and references (genome, annotation)
  - [ ] check chromosomes
  - [ ] check mate ids
  - [ ] add trimming and quality filtering step (e.g. cutadapt)
- [ ] allow for custom arguments in tool templates
- [ ] more flexible tool configuration
- [ ] add `docker pull` command as `beforeScript` to config
- [x] add pipeline to CI server (circle)
- [ ] integrate Kallisto, Sailfish and maybe other mappers
- [ ] tune mapping steps for different quantification modes: only produce required bams
- [ ] autodetect sjdbOverHang for STAR index
- [ ] generate gff files from RSEM quantification
- [ ] use `ext` scope for task properties