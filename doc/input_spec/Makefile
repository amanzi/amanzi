#
#  Makefile: generating HTML and PDF from RST
#
#  NB: not a good solution, just an initial hack.

CBpath =./code-block
PY_CSS = --stylesheet=$(CBpath)/pygments-default.css
PY_STY = --stylesheet=$(CBpath)/pygments-default.sty

SPHINXOPTS  =
SPHINXBUILD = sphinx-build
SOURCEDIR   = ./
BUILDDIR    = build

SRC = \
	AmanziInputSpec-u.rst \
	AmanziInputSpec-s.rst \
        AmanziNativeSpec_dev.rst.pre \
        AmanziNativeSpec_dev.rst

RST = $(SRC:%.rst.pre=%.rst)

HTML = $(SRC:%.rst=%.html)

PDF = $(SRC:%.rst=%.pdf)

# generic rules for targets

%.html : %.rst
	$(CBpath)/rst2html-pygments $(PY_CSS) $< $@ 

%.tex : %.rst
	$(CBpath)/rst2latex-pygments $(PY_STY) $< $@

%.pdf : %.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<

%.rst:
	@./AmanziNativeSpec.py $< $@

sphinx:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

rst : $(RST)

html : $(HTML)

all : rst sphinx

pdf : $(PDF)

clean-rst : 
	-/bin/rm -f AmanziNativeSpec_dev.rst

clean-html : 
	-/bin/rm -f *.html

clean-pdf :
	-/bin/rm -f *.tex *.pdf
	-/bin/rm -f *.aux *.toc *.log *.out

clean-sphinx : 
	-/bin/rm -rf build

clean : clean-html clean-pdf clean-rst clean-sphinx
	-/bin/rm -f *~

distclean : clean

