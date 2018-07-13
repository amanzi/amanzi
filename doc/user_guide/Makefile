# Makefile for Sphinx documentation
#

AMANZI_UG_DIR = $(CURDIR)
AMANZI_PYTHON = $(subst doc/user_guide,tools/amanzi_xml,$(CURDIR)):$(subst doc/user_guide,tools/prettytable,$(CURDIR))

ifndef PYTHONPATH
  NEW_PYTHONPATH=$(AMANZI_PYTHON)
else
  NEW_PYTHONPATH=$(AMANZI_PYTHON):$(PYTHONPATH)
endif

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = export PYTHONPATH=$(NEW_PYTHONPATH); sphinx-build
PAPER         =
BUILDDIR      = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .


.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  html-ssl   to make HTML that use SSL link for MathJax"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"
	@echo ""
	@echo "Environment Check "
	@echo "  AMANZI_INSTALL_DIR = $(AMANZI_INSTALL_DIR)"
	@echo "  AMANZI_MPI_EXEC    = $(AMANZI_MPI_EXEC)"
	@echo "  AMANZI_MPI_NP      = $(AMANZI_MPI_NP)"
	@echo "  AMANZI_UG_DIR      = $(AMANZI_UG_DIR)"
	@echo "  AMANZI_PYTHON      = $(AMANZI_PYTHON)"
	@echo "  PYTHONPATH         = $(PYTHONPATH)"
	@echo "  NEW_PYTHONPATH     = $(NEW_PYTHONPATH)"
	@echo ""
	@echo "Use the variable AMANZI_SECTION_OPTIONS to include optional sections of the"
	@echo "User Guide in your build,"
	@echo ""
	@echo "   AMANZI_SECTION_OPTIONS = $(AMANZI_SECTION_OPTIONS)"
	@echo ""
	@echo "This variables are passed to the collect.py script ... "
	@python collect.py --help
	@echo ""
	@echo "Use the variable AMANZI_INCLUDE_TODOS=(True|False) to include TODO elements,"
	@echo ""
	@echo "   AMANZI_INCLUDE_TODOS = $(AMANZI_INCLUDE_TODOS)"
	@echo ""
	@echo "Both of these variables are read by the Sphinx conf.py script."

clean:
	@echo "Removing the build directory: $(BUILDDIR)" ;\
	rm -rf $(BUILDDIR)/*  ;\
	echo "Removing the section subdirectories: verification, tutorial, benchmarking, install, mycase" ;\
	rm -rf verification tutorial benchmarking install mycase  ;\
	echo "Removing the generated ./index.rst"   ;\
	rm -f ./index.rst

run_tests:
	@echo "Setting up guide indexes, and running any tests for this build ..."
	export PYTHONPATH=$(NEW_PYTHONPATH);\
	python collect.py $(AMANZI_SECTION_OPTIONS) --run-tests
	@echo "TESTs finished."
	@echo

html:   run_tests
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

html-norun:
	export PYTHONPATH=$(NEW_PYTHONPATH); python collect.py $(AMANZI_SECTION_OPTIONS)
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

html-ssl:
	make MATHJAX_SSL=1 html

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(BUILDDIR)/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/Amanzi.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/Amanzi.qhc"

devhelp:
	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(BUILDDIR)/devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/Amanzi"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/Amanzi"
	@echo "# devhelp"

epub:
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf:
	python collect.py $(AMANZI_SECTION_OPTIONS) 
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

text:
	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man:
	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(BUILDDIR)/man
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext:
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."
