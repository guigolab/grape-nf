.PHONY: changelog docs test clean-docs clean

changelog:
	@git log --pretty=%s $$(git describe --abbrev=0)..HEAD

docs: docs/Makefile
	@echo -n "== $$(tput setaf 3; tput bold)Updating GRAPE documentation...$$(tput sgr0)"
	@$(MAKE) -C docs deploy
	@echo "$$(tput setaf 2; tput bold)done$$(tput sgr0)"

test: test.sh
	@PROFILE=testStarRsem ./test.sh

clean-docs: docs/Makefile
	@$(MAKE) -C docs clean
	
clean:
	@rm -rf work nextflow .nextflow* checksum pipeline.db
	