EPStoPDF = epstopdf

ALL: paper.ps paper.pdf 

eps_source = $(wildcard plots/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

paper.dvi: paper.tex ../refs.bib $(eps_source)
	latex paper.tex < /dev/null
	bibtex paper < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null

paper.pdf: paper.tex ../refs.bib $(pdf_source)
	pdflatex paper.tex < /dev/null
	bibtex paper < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null

pdf:	paper.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) *.dvi 
	$(RM) *.blg *.log *.bib
	$(RM) *.aux *.ps *.bbl
	$(RM) paper.pdf
	$(RM) *~

.PHONY: clean
