.PHONY: all html check_deploy deploy clean deepclean

GEMS = env
HTML_DIR = _html
PDF_DIR = _pdf
PHP_DIR = _php
DOCS = quickstart
INDEX = quickstart
PHP_HEADER = header.php
DEPLOY_SRC = _deploy
ifndef DEPLOY_EXT
	DEPLOY_EXT = php
endif

all: html

html: setup $(DOCS:=.html)
pdf: setup $(DOCS:=.pdf)

php: html $(DOCS:=.php)

$(DOCS:=.html): %.html: src/%.adoc
	@GEM_HOME=$(GEMS) $(GEMS)/bin/asciidoctor $^ -o $(HTML_DIR)/$@

$(DOCS:=.pdf): %.pdf: src/%.adoc
	@GEM_HOME=$(GEMS) $(GEMS)/bin/asciidoctor-pdf $^ -o $(PDF_DIR)/$@

$(DOCS:=.php): %.php: $(HTML_DIR)/%.html
	@mkdir -p $(PHP_DIR)
	@cat src/$(PHP_HEADER) $^ >> $(PHP_DIR)/$@

setup: $(GEMS)/bin/asciidoctor $(GEMS)/bin/asciidoctor-pdf

$(GEMS)/bin/asciidoctor: $(GEMS)
	@GEM_HOME=$(GEMS) $(GEMS)/bin/bundle --path=$(GEMS) --binstubs=$(GEMS)/bin

$(GEMS): $(GEMS)/bin/bundle
$(GEMS)/bin/bundle:
	@GEM_HOME=$(GEMS) gem install bundler

prepare_deploy: $(DEPLOY_EXT)
	@mkdir -p $(DEPLOY_SRC)
	@ln -f src/password.php $(DEPLOY_SRC)
	@for doc in $(DOCS:=.$(DEPLOY_EXT)); do \
		[[ $$doc == $(INDEX).$(DEPLOY_EXT) ]] && DEST=index.$(DEPLOY_EXT) || DEST=$$doc; \
		ln -f $(PHP_DIR)/$$doc $(DEPLOY_SRC)/$$DEST; \
	done

check_deploy:
ifndef DEPLOY_DEST
	$(error Undefined variable DEPLOY_DEST)
endif

deploy: prepare_deploy check_deploy
	rsync -a --delete --exclude=".*" $(DEPLOY_SRC)/ $(DEPLOY_DEST)/

serve:
	php -S 127.0.0.1:8000 -t .

clean:
	rm -rf $(HTML_DIR) $(PDF_DIR) $(DEPLOY_SRC) $(PHP_DIR)

deepclean: clean
	[[ $(GEMS) == 'env' ]] && rm -rf $(GEMS) 
	rm -rf bin .bundle
