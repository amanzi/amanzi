#
#  Makefile: generating HTML and PDF from RST
#
#  NB: not a good solution, just an initial hack.

CBpath =./code-block/
PY_CSS = --stylesheet=$(CBpath)/pygments-default.css
PY_STY = --stylesheet=$(CBpath)/pygments-default.sty

html: AmanziInputSpec.rst
	$(CBpath)/rst2html-pygments $(PY_CSS) AmanziInputSpec.rst AmanziInputSpec.html 

pdf: AmanziInputSpec.rst
	$(CBpath)/rst2latex-pygments $(PY_STY) AmanziInputSpec.rst AmanziInputSpec.tex; \
	pdflatex AmanziInputSpec.tex; \
	pdflatex AmanziInputSpec.tex; \
	pdflatex AmanziInputSpec.tex

native-html: AmanziNativeSpecV2.rst
	$(CBpath)/rst2html-pygments $(PY_CSS) AmanziNativeSpecV2.rst AmanziNativeSpecV2.html


all: html pdf

clean-html: 
	/bin/rm -f AmanziInputSpec.html

clean-pdf:
	/bin/rm -f AmanziInputSpec.tex AmanziInputSpec.pdf 
	/bin/rm -f *.aux *.toc *.log *.out

clean: clean-html clean-pdf

distclean: clean
