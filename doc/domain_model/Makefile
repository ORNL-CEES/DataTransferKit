TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: 
	$(TEX) domain_model
	$(BIB) domain_model
	$(TEX) domain_model
	$(TEX) domain_model
	$(PDF) domain_model

clean:
	-rm -f *~
	-rm -f *.log
	-rm -f *.aux
	-rm -f *.dvi
	-rm -f *.rel
	-rm -f *.bbl
	-rm -f *.blg
	-rm -f *.out
	-rm -f *.toc *.lof *.lot
	-rm -f *.bak
	-rm -f TAGS

realclean: clean
	-rm -f domain_model.pdf
