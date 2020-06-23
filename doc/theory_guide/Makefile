
SHELL=/bin/sh

FILE=theory
TEXFILE=$(FILE).tex
PDFFILE=$(FILE).pdf
#BIBFILE=$(FILE).bib

#
# texfile dependencies
#
TEXFILES=$(shell find . -name "*.tex" | sed -e 's/\.\///' )
#BIBFILES=$(shell find ./bib -name "*\.bib" | sort -n )
AUXFILES=$(shell ls bu*.aux)

$(PDFFILE): $(TEXFILES) $(BIBFILE)
	TEXINPUTS=.//:;\
	pdflatex $(TEXFILE);\
	bibtex $(FILE);\
	pdflatex $(TEXFILE);\
	pdflatex $(TEXFILE)


$(BIBFILE): $(BIBFILES)
	/bin/rm -f $(BIBFILE); cat $(BIBFILES) > $(BIBFILE)

first:  $(TEXFILES)
	TEXINPUTS=.//:;\
	pdflatex $(TEXFILE)

clean:
	/bin/rm -f $(PDFFILE)

distclean:
	/bin/rm -f $(PDFFILE) $(BIBFILE)
	/bin/rm -f  *.aux *.dvi *.log *.bbl *.blg *.toc *.lof *.lot *.out
