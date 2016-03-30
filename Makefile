docs: docs/Makefile
	@echo -n "== $$(tput setaf 3; tput bold)Updating GRAPE documentation...$$(tput sgr0)"
	@$(MAKE) -C docs deploy
	@echo "$$(tput setaf 2; tput bold)done$$(tput sgr0)"

test: test.sh
	@./test.sh

clean-docs: docs/Makefile
	@$(MAKE) -C docs clean
	
clean-test:
	@rm -rf work nextflow .nextflow* checksum pipeline.db

clean-all: clean-docs clean-test
