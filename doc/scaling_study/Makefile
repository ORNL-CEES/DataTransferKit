TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: 
	$(PDFTEX) dtk_scaling_study
	$(BIB) dtk_scaling_study
	$(PDFTEX) dtk_scaling_study
	$(PDFTEX) dtk_scaling_study

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
	-rm -f dtk_scaling_study.pdf
