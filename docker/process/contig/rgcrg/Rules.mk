RGCRG_VER = 0.1

.contig_rgcrg: bedtools sambamba
.contig_rgcrg: RGCRG_ARGS = \
		BEDTOOLS_VER=$(BEDTOOLS_VER) \
		SAMBAMBA_VER=$(SAMBAMBA_VER)