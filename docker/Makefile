-include common.mk
-include $(PROCS:%=process/%/Rules.mk)

VERBOSE ?= 1
PROC_TOOLS := $(foreach p,$(PROCS),$(foreach t,$($(shell echo $(p) | tr a-z A-Z)_TOOLS),.$(p)_$(t)))

.PHONY: processes tools $(TOOLS) $(PROC_TOOLS)

processes: $(PROC_TOOLS)

tools: $(TOOLS)

logs/tools:
	@mkdir -p logs/tools

$(TOOLS): TOOL=$(shell echo $@ | cut -d- -f1| tr a-z A-Z)
$(TOOLS): logs/tools
	@[ $(VERBOSE) ] && echo Building container for $@ $($(TOOL)_VER) || true
	@echo docker build -t grapenf/$@:$($(TOOL)_VER) tool/$@ #> logs/$@.log

$(PROC_TOOLS): INFO=$(subst _, ,$@)
$(PROC_TOOLS): PROC=$(word 1,$(INFO))
$(PROC_TOOLS): PROC_TOOL=$(word 2,$(INFO))
$(PROC_TOOLS): TOOL=$(shell echo $(PROC_TOOL) | tr a-z A-Z)
$(PROC_TOOLS):
	@[ $(VERBOSE) ] && echo Building container for $@ $($(TOOL)_VER) || true
	@echo docker build $($(TOOL)_ARGS:%=--build-arg %) \
		-t grapenf/$(PROC):$(PROC_TOOL)-$($(TOOL)_VER) process/$(PROC)/$(PROC_TOOL)