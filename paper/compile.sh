# alternatively -t CFONT
latexdiff -t BOLD -s COLOR --no-del paper_submitted.tex paper.tex > paper_diff.tex
latexmk -f
