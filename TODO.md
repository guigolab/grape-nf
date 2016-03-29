# Things to do

- [ ] add preprocessing steps for input (fastq) and references (genome, annotation)
  - [ ] check chromosomes
  - [ ] check mate ids
- [ ] allow for custom arguments in tool templates
- [ ] more flexible tool configuration
- [ ] add `docker pull` command as `beforeScript` to config
- [ ] add pipeline to CI server (circle)
- [ ] integrate Kallisto, Sailfish and maybe other mappers
- [ ] tune mapping steps for different quantification modes: only produce required bams
- [ ] autodetect sjdbOverHang for STAR index
- [ ] generate gff files from RSEM quantification
