import panflute as pf

def equation_environment(math_elem, doc):
    # Check for paragraphs containing only a DisplayMath element
    if isinstance(math_elem, pf.Math) and math_elem.format == "DisplayMath":
        # Replace with raw LaTeX equation environment
        tex = f"\\begin{{equation}}{math_elem.text}\\end{{equation}}"
        return pf.RawInline(tex, format="latex")


def main(doc=None):
    return pf.run_filter(equation_environment, doc=doc)

if __name__ == "__main__":
    main()
