MAPPING_TOOLS = gem star
MAPPING_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: $(MAPPING_TOOLS)

-include $(MAPPING_TOOLS:%=$(MAPPING_DIR)/%/Rules.mk)