set -e 

pandocoptions="--defaults=pandoc_params.yaml"

pandoc $pandocoptions -o abstract.tex abstract.md
pandoc $pandocoptions -o introduction.tex introduction.md
pandoc $pandocoptions -o data.tex data.md
pandoc $pandocoptions -o methods.tex methods.md
pandoc $pandocoptions -o results.tex results.md
pandoc $pandocoptions -o discussion.tex discussion.md
pandoc $pandocoptions -o acknowledgements.tex acknowledgements.md
pandoc $pandocoptions -o research_acknowledgements.tex research_acknowledgements.md

pandoc $pandocoptions -o extra_density.tex extra_density.md --top-level-division=chapter
pandoc $pandocoptions -o extra_rv_models.tex extra_rv_models.md --top-level-division=chapter
pandoc $pandocoptions -o extra_numerical_convergence.tex extra_numerical_convergence.md --top-level-division=chapter
pandoc $pandocoptions -o extra_results.tex extra_results.md --top-level-division=chapter

