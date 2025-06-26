set -e 
rsync ../dwarfs_figures/thesis/plots/figures/ figures -av

# pandocoptions="--filter ./table_shortcaption.py --filter pandoc-crossref --citeproc -M --top-level-division=section --bibliography=main.bib --natbib -M --autoSectionLabels=false -M --autoEqnLabels=true --filter ./table.py --filter ./md_figure.py --number-sections --filter equation.py"
pandocoptions="--defaults=pandoc_params.yaml"

pandoc $pandocoptions -o abstract.tex abstract.md
pandoc $pandocoptions -o introduction.tex introduction.md
pandoc $pandocoptions -o data.tex data.md
pandoc $pandocoptions -o methods.tex methods.md
pandoc $pandocoptions -o sculptor.tex sculptor.md

# pandoc $pandocoptions -o ursa_minor.tex ursa_minor.md
# pandoc $pandocoptions -o conclusion.tex conclusion.md

# pandoc $pandocoptions -o appendix.tex appendix.md --top-level-division=chapter

