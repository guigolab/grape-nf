ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

TOOLS := bamstats bedtools flux-capacitor gemtools kentutils rsem rseqc sambamba samtools star
PROCS = contig quantification

GLIBC_VER := 2.25

BAMSTATS_VER := 0.3.2
BEDTOOLS_VER := 2.19.1
FLUX_VER := 1.6.1
GEMTOOLS_VER := 1.7.1
KENTUTILS_VER := 308
RSEM_VER := 1.2.21
RSEQC_VER := 2.6.4
SAMBAMBA_VER := 0.7.0
SAMTOOLS_VER := 1.3.1
STAR_VER := 2.4.0j

RGCRG_VER := 0.1