FLAGS = -shell-escape 

all: myrefs.bib SLGF_Journal.tex
	pdflatex $(FLAGS)  SLGF_Journal.tex 
	bibtex SLGF_Journal.aux 
	pdflatex $(FLAGS) SLGF_Journal.tex 
	pdflatex $(FLAGS) SLGF_Journal.tex

differences: myrefs.bib SLGF_Journal.tex
	latexdiff-vc --svn -r SLGF_Journal.tex
	pdflatex $(FLAGS)  SLGF_Journal-diff.tex 
	bibtex SLGF_Journal-diff.aux 
	pdflatex $(FLAGS) SLGF_Journal-diff.tex 
	pdflatex $(FLAGS) SLGF_Journal-diff.tex

