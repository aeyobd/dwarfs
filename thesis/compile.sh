set -e 
rsync ../dwarfs_figures/thesis/plots/figures/ figures -av

pandocoptions="--filter pandoc-crossref --citeproc -M --top-level-division=section --bibliography=main.bib --natbib -M --autoSectionLabels=false -M --autoEqnLabels=true --filter ./table.py --filter ./md_figure.py"

pandoc $pandocoptions -o abstract.tex abstract.md
pandoc $pandocoptions -o introduction.tex introduction.md
pandoc $pandocoptions -o data.tex data.md
pandoc $pandocoptions -o methods.tex methods.md
pandoc $pandocoptions -o numerical_convergence.tex numerical_convergence.md


latexmk -f
