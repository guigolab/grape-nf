ifndef DEPLOYDIR
	DEPLOYDIR = /users/rg/epalumbo/public_html/grape/docs
endif

docs: $(DEPLOYDIR)
	@echo -n "== $$(tput setaf 3; tput bold)Updating http://genome.crg.es/~epalumbo/grape...$$(tput sgr0)"
	@source /software/rg/el6.3/virtualenvs/indexfile.devel/bin/activate && make -C docs deploy DEPLOYDIR=$(DEPLOYDIR) &> .html.log 
	@echo "$$(tput setaf 2; tput bold)done$$(tput sgr0)"


test: compare.sh
	@./compare.sh
