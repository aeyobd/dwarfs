# alternatively -t CFONT
latexdiff -t BOLD --no-del paper_rev1.tex paper.tex > paper_diff.tex
latexmk -f
