#
#  Makefile: generating HTML and PDF from RST
#
#  NB: not a good solution, just an initial hack.

CBpath =./code-block
PY_CSS = --stylesheet=$(CBpath)/pygments-default.css
PY_STY = --stylesheet=$(CBpath)/pygments-default.sty

SRC = \
	ATSNativeSpec_dev.rst.in

RST = $(SRC:%.rst.in=%.rst)

HTML = $(SRC:%.rst.in=%.html)

PDF = $(SRC:%.rst.in=%.pdf)

# generic rules for targets

%.rst:
	@./generate_ATSNativeSpec.py $< $@

%.html : %.rst
	$(CBpath)/rst2html-pygments $(PY_CSS) $< $@ 

%.tex : %.rst
	$(CBpath)/rst2latex-pygments $(PY_STY) $< $@

%.pdf : %.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<

all : rst html pdf

rst : $(RST)

html : $(HTML)

pdf : $(PDF)

clean-rst :
	-/bin/rm -f $(RST)

clean-html : 
	-/bin/rm -f *.html

clean-pdf :
	-/bin/rm -f *.tex *.pdf
	-/bin/rm -f *.aux *.toc *.log *.out

clean : clean-rst clean-html clean-pdf
	-/bin/rm -f *~

distclean : clean
