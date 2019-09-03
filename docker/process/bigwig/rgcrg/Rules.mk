RGCRG_VER = 0.1

.bigwig_rgcrg: bedtools kentutils sambamba
.bigwig_rgcrg: RGCRG_ARGS = \
		BEDTOOLS_VER=$(BEDTOOLS_VER) \
		KENTUTILS_VER=$(KENTUTILS_VER) \
		SAMBAMBA_VER=$(SAMBAMBA_VER)