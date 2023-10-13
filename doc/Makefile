#Makefile for generating the documentation.

filename1=SFT1D

PDF_FILES=usermanual.pdf 

pdf:
	pdflatex ${filename1}.tex > /dev/null
	bibtex ${filename1} > /dev/null
	#makeindex ${filename1} > /dev/null
	pdflatex ${filename1}.tex > /dev/null
	pdflatex ${filename1}.tex > /dev/null
	mv ${filename1}.pdf ${PDF_FILES}

clean:
	rm -f *.out *.aux *.log *.toc *.ind *.idx *.ilg *.dvi *.ps *~ *.bbl *.blg

