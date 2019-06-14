INFEREXP_TOOLS = rseqc
INFEREXP_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: $(INFEREXP_TOOLS)

-include $(INFEREXP_TOOLS:%=$(INFEREXP_DIR)/%/Rules.mk)