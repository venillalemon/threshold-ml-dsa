LATEXMK ?= latexmk
DOCUMENT := slot-protocol-figure_3

.PHONY: all clean distclean
all: $(DOCUMENT).pdf

$(DOCUMENT).pdf: $(DOCUMENT).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -file-line-error $<

clean:
	$(LATEXMK) -c $(DOCUMENT).tex

distclean:
	$(LATEXMK) -C $(DOCUMENT).tex
