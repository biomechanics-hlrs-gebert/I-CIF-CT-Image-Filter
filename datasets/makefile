# --------------------------------------------------------------------------------------------------
# Makefile to postprocess CTIF data sets.
#
# Author:    Johannes Gebert, gebert@hlrs.de, HLRS, NUM
# Date:      25.12.2021
# Last edit: 25.12.2021
#
# For use of make visit: https://www.gnu.org/software/make/
# --------------------------------------------------------------------------------------------------
LATEX_CMD=lualatex
BIBER_CMD=biber
#
auxsuffixes=*.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl *.bcf *.nav *.run.xml *.synctex.gz *.snm *.fls *.fdb_latexmk
#
all: ${tex}
ifeq ($(tex),)
	@echo "-----------------------------------------------"
	@echo "-- Usage:"
	@echo "--      make tex=<*.tex file>"
	@echo "-----------------------------------------------"
endif
#
pdf: ${tex}
#
final: ${tex}
	rm -f ${auxsuffixes}
#
clean:
	rm -f ${auxsuffixes}
#
${tex}: clean ${tex}
	$(LATEX_CMD) ${tex}