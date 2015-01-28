EPStoPDF = epstopdf

guide = userguide.tex

ALL: userguide.ps userguide.pdf 

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

userguide.dvi: userguide.tex $(eps_source)
	latex userguide.tex < /dev/null
	bibtex userguide
	latex userguide.tex < /dev/null
	latex userguide.tex < /dev/null
	latex userguide.tex < /dev/null

userguide.pdf: userguide.tex $(pdf_source)
	pdflatex userguide.tex < /dev/null
	bibtex userguide
	pdflatex userguide.tex < /dev/null
	pdflatex userguide.tex < /dev/null
	pdflatex userguide.tex < /dev/null

pdf:	userguide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) userguide.dvi
	$(RM) userguide.blg userguide.log userguide.toc userguide.out
	$(RM) userguide.aux userguide.ps userguide.bbl userguide.pdf
	$(RM) *~

.PHONY: clean
