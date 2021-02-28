all:	WFN1C.pdf
# WFN1.pdf WFN1_supplementary.pdf

WFN1C.pdf: WFN1C.tex WFN1C.bib WFN1C_alt.bib
	pdflatex WFN1C.tex
	bibtex WFN1C
	pdflatex WFN1C.tex
	pdflatex WFN1C.tex

WFN1.pdf: WFN1.tex WFN1.bib WFN1_alt.bib
	pdflatex WFN1.tex
	bibtex WFN1
	pdflatex WFN1.tex
	pdflatex WFN1.tex

WFN1_supplementary.pdf: WFN1_supplementary.tex WFN1.bib WFN1_alt.bib
	pdflatex WFN1_supplementary.tex
	bibtex WFN1_supplementary
	pdflatex WFN1_supplementary.tex
	pdflatex WFN1_supplementary.tex

