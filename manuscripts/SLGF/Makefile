FLAGS = -shell-escape 

all: myrefs.bib SLGF_RECOMB.tex
	pdflatex $(FLAGS)  SLGF_RECOMB.tex 
	bibtex SLGF_RECOMB.aux 
	pdflatex $(FLAGS) SLGF_RECOMB.tex 
	pdflatex $(FLAGS) SLGF_RECOMB.tex

differences: myrefs.bib SLGF_RECOMB.tex
	latexdiff-vc --svn -r SLGF_RECOMB.tex
	pdflatex $(FLAGS)  SLGF_RECOMB-diff.tex 
	bibtex SLGF_RECOMB-diff.aux 
	pdflatex $(FLAGS) SLGF_RECOMB-diff.tex 
	pdflatex $(FLAGS) SLGF_RECOMB-diff.tex

abstract: SLGF_RECOMB_ABSTRACT.tex
	pdflatex $(FLAGS)  SLGF_RECOMB_ABSTRACT.tex
	bibtex SLGF_RECOMB_ABSTRACT.aux
	pdflatex $(FLAGS) SLGF_RECOMB_ABSTRACT.tex
	pdflatex $(FLAGS) SLGF_RECOMB_ABSTRACT.tex
