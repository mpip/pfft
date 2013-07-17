
# add all source files here
manual_FILES = \
	manual.tex \
	intro.tex \
	install.tex \
	fortran.tex \
	shortcuts.tex \
	tutorial.tex \
	features.tex \
	bibliography.bib \
	code/manual_min_c2c.c

EXTRA_DIST = $(manual_FILES)
if HAVE_LATEX
EXTRA_DIST += manual.pdf
endif

.PHONY: doc manual

doc: manual
manual: manual.pdf
manual.pdf: $(manual_FILES) latex-manual.sh latexinc.sh
if HAVE_LATEX
	sh latex-manual.sh
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the User's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

CLEANFILES = manual.pdf

MOSTLYCLEANFILES = \
	*.aux *.aux.bak\
	manual.bbl manual.blg \
	manual.idx manual.idx.bak manual.ilg manual.ind\
	manual.log manual.out manual.toc manual.lof manual.lot manual.tdo
