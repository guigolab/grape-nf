.PHONY: all html check_deploy deploy clean deepclean

GEMS = env
DEPLOY_LIST = deploy-list.txt
HTML_DIR = _html
PDF_DIR = _pdf
DOCS = quickstart workflow
SECURE_HEADER = header.php

all: html

html: setup $(DOCS:=.html)
pdf: setup $(DOCS:=.pdf)

$(DOCS:=.html): %.html: src/%.adoc
	@GEM_HOME=$(GEMS) $(GEMS)/bin/asciidoctor $^ -o $(HTML_DIR)/$@

$(DOCS:=.pdf): %.pdf: src/%.adoc
	@GEM_HOME=$(GEMS) $(GEMS)/bin/asciidoctor-pdf $^ -o $(PDF_DIR)/$@

setup: $(GEMS)/bin/asciidoctor $(GEMS)/bin/asciidoctor-pdf

$(GEMS)/bin/asciidoctor: $(GEMS)
	@GEM_HOME=$(GEMS) $(GEMS)/bin/bundle --path=$(GEMS) --binstubs=$(GEMS)/bin

$(GEMS): $(GEMS)/bin/bundle
$(GEMS)/bin/bundle:
	@GEM_HOME=$(GEMS) gem install bundler

secure: html $(DOCS:=.php)

$(DOCS:=.php): %.php: $(HTML_DIR)/%.html
	@cat $(SECURE_HEADER) $^ >> $@

$(DEPLOY_LIST):
	@for doc in $(DOCS:=.php); do \
		echo $$doc >> $@; \
	done
	@echo '.htpasswd' >> $@

check_deploy:
ifndef DEPLOY_DIR
	$(error Undefined variable DEPLOY_DIR)
endif

deploy: html $(DEPLOY_LIST) check_deploy
	rsync -a --files-from=$(DEPLOY_LIST) . $(DEPLOY_DIR)

serve:
	php -S 127.0.0.1:8000 -t .

clean:
	rm -rf $(HTML_DIR) $(PDF_DIR) $(DEPLOY_LIST) $(DOCS:=.php)

deepclean: clean
	rm -rf $(GEMS) bin .bundle
