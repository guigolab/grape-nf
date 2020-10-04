.quantification_rsem: rsem sambamba
.quantification_rsem: TAGS = $(NS)/txindex:rsem-${RSEM_VER}
.quantification_rsem: RSEM_ARGS = \
		RSEM_VER=$(RSEM_VER) \
		SAMBAMBA_VER=$(SAMBAMBA_VER)