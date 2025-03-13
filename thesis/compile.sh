set -e 
rsync ../dwarfs_figures/thesis/plots/figures/ figures -av

pandocoptions="--filter pandoc-crossref --citeproc -M --top-level-division=section --bibliography=main.bib --natbib"

pandoc $pandocoptions -o introduction.tex introduction.md
latexmk -xelatex -f main.tex
