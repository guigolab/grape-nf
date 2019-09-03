CONTIG_TOOLS = rgcrg
CONTIG_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: $(CONTIG_TOOLS)

-include $(CONTIG_TOOLS:%=$(CONTIG_DIR)/%/Rules.mk)