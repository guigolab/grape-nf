BIGWIG_TOOLS = rgcrg star
BIGWIG_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: $(BIGWIG_TOOLS)

-include $(BIGWIG_TOOLS:%=$(BIGWIG_DIR)/%/Rules.mk)