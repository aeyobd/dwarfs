# alternatively -t CFONT
latexdiff -t BOLD --no-del paper_submitted.tex paper.tex > paper_diff.tex
latexmk -f
