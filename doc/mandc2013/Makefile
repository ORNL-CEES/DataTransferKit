TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: stencil.pdf
	$(PDFTEX) slattery_wilson_pawlowski
	$(BIB) slattery_wilson_pawlowski
	$(PDFTEX) slattery_wilson_pawlowski
	$(PDFTEX) slattery_wilson_pawlowski
	$(PDFTEX) slattery_wilson_pawlowski

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
	-rm -f *.pdftex_t
	-rm -f TAGS
	-rm -f *.pdf
