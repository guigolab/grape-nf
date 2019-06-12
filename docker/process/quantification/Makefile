RSEM_VER = 1.2.21
SAMBAMBA_VER = 0.7.0
FLUX_VER = 1.6.1
SAMTOOLS_VER = 1.3.1

TOOLS = rsem flux

.PHONY: $(TOOLS)

default: $(TOOLS)

rsem: RSEM_ARGS = \
		RSEM_VER=$(RSEM_VER) \
		SAMBAMBA_VER=$(SAMBAMBA_VER)

flux: FLUX_ARGS = \
		FLUX_VER=$(FLUX_VER) \
		SAMTOOLS_VER=$(SAMTOOLS_VER)

$(TOOLS): TOOL=$(shell echo $@ | tr a-z A-Z)

$(TOOLS):
	@docker build $($(TOOL)_ARGS:%=--build-arg %) \
		-t grapenf/quantification:$@-$($(TOOL)_VER) $@