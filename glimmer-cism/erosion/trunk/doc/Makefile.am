if DO_DOC_PDF
 pdf_doc = erosion.pdf
endif
if DO_DOC_PS
 ps_doc = erosion.ps
endif
doc_DATA = $(ps_doc) $(pdf_doc)
docdir = $(prefix)/doc
dist_doc_DATA = $(doc_DATA)
SUBDIRS = erosion common

EXTRA_DIST = erosion.tex  erosion_pdf.tex  erosion_ps.tex erosion.bib

CLEANFILES = *~ *.aux *.dvi *.log *.out *.toc

erosion.pdf:	erosion_pdf.dvi
		$(DVIPDFT) -o erosion.pdf erosion_pdf.dvi

erosion.ps:	erosion_ps.dvi
		dvips -o erosion.ps erosion_ps.dvi

.tex.dvi::
	latex $<
	bibtex $*
	latex $<
