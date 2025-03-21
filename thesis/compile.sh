set -e 
rsync ../dwarfs_figures/thesis/plots/figures/ figures -av

pandocoptions="--filter pandoc-crossref --citeproc -M --top-level-division=section --bibliography=main.bib --natbib --filter ./table.py --filter ./md_figure.py"

pandoc $pandocoptions -o abstract.tex abstract.md
pandoc $pandocoptions -o introduction.tex introduction.md
latexmk -xelatex -f main.tex
