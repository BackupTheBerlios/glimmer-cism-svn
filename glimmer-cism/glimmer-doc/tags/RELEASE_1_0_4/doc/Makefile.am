if DO_DOC_PDF
 pdf_doc = glimmer.pdf
endif
if DO_DOC_PS
 ps_doc = glimmer.ps
endif
if HAVE_PYCF
 pycf_dir = tg
endif
doc_DATA = $(ps_doc) $(pdf_doc)
docdir = $(prefix)/doc
dist_doc_DATA = $(doc_DATA)
SUBDIRS = dg ug num tut common $(pycf_dir)

EXTRA_DIST = dg.tex dg_pdf.tex dg_ps.tex glimmer.tex glimmer_pdf.tex \
             glimmer_ps.tex num.tex num_pdf.tex num_ps.tex ug.tex ug_pdf.tex \
             ug_ps.tex tg_pdf.tex tg_ps.tex tg.tex tut_pdf.tex tut_ps.tex tut.tex glimmer.bib

CLEANFILES = *~ *.aux *.dvi *.log *.out *.toc *.bbl *blg *.html \
             glimmer_html.4ct  glimmer_html.4tc  glimmer_html.css  glimmer_html.idv  \
             glimmer_html.lg  glimmer_html.tmp  glimmer_html.xref

code_doc.ps:	code_doc.dvi
		dvips -o code_doc.ps code_doc.dvi

code_doc.tex:	../src/fortran/*.f90
		f90docp.py -t "GLIMMER Code Documentation" -a "Glimmer Developer Team" $^ > $@
www::	glimmer.tex
	install -d www
	htlatex glimmer_html "html,3,info" "" " -d./www/ -m 644 "
	$(BIBTEX) glimmer_html
	htlatex glimmer_html "html,3,info" "" " -d./www/ -m 644 "
	find $(SUBDIRS) -name '*.png' -exec install -m 644 -D {} ./www/{} \;

.tex.dvi::
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<
	$(LATEX) $<

%.ps:	%_ps.dvi
	dvips -o $@ $<

%.pdf:	%_pdf.dvi
	$(DVIPDFT) -o $@ $<
