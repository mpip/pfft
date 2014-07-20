
# add all source files here
manual_FILES = \
	manual.tex \
	intro.tex \
	install.tex \
	fortran.tex \
	shortcuts.tex \
	tutorial.tex \
	features.tex \
	bibliography.bib

EXTRA_DIST = $(manual_FILES)
CLEANFILES =
MOSTLYCLEANFILES =

.PHONY: doc manual

doc: manual

if ENABLE_DOC
if HAVE_LATEX
manual: manual.pdf
else
manual:
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the User's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif
else
manual:
	@echo "Skip build of documentation because of --disable-doc."
endif

if ENABLE_DOC
if HAVE_LATEX
%.tex: ../tests/%.ctex
	cat $< \
	    | sed '/^\\del§[^§]*§$$/d;         				    # erase whole line if \del spans from start to end' \
	    | sed '/\\del§/{s,\\del§[^§]*§,,g} 				    # apply single-line \del' \
	    | sed '/^\\del§$$/,/^§$$/d         				    # apply multi-line \del' \
	    | sed ':a s,\\red§\([^§]*[^\\]\)_\([^§]*\)§,\\red§\1\\_\2§,; ta # replace _ by \_ in single-line \red' \
	    | sed '/^\\red§$$/,/^§$$/{s/_/\\_/g}			    # replace _ by \_ in multi-line \red' \
	    | sed 's/\\red§/(\*@\\color{red}/g; s/§/@\*)/g		    # convert \red tags into LaTeX (single- and multi-line)' \
	    > $(builddir)/$@

# Following files are not part of EXTRA_DIST, since they are built from .ctex files
gentex_FILES = \
	man_fftw_3don1d.tex

manual.pdf: $(manual_FILES) $(gentex_FILES) latex-manual.sh latexinc.sh
	sh latex-manual.sh

EXTRA_DIST  += manual.pdf 
CLEANFILES  += manual.pdf
MOSTLYCLEANFILES += \
	*.aux *.aux.bak\
	manual.bbl manual.blg \
	manual.idx manual.idx.bak manual.ilg manual.ind\
	manual.log manual.out manual.toc manual.lof manual.lot manual.tdo \
	$(gentex_FILES)
endif
endif

