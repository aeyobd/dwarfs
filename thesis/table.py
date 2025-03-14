# from discussion https://gist.github.com/wgroeneveld/9dbeb0d0b60c6cb5d8dfe9b938c5e94e
import sys

from panflute import *

def replace_longtables_with_tabular(elem, doc):
    try:
        def get_text(item):
            if isinstance(item, Str):
                return item.text.replace('_', r'\_')
            elif isinstance(item, Plain):
                return ''.join([get_text(i) for i in item.content])
            elif isinstance(item, ListContainer):
                return ''.join([get_text(i) for i in item])
            elif isinstance(item, Code):
                return '\\verb|' + (item.text) + '|'
            elif isinstance(item, Math):
                return '$' + (item.text) + '$'
            elif isinstance(item, RawInline):
                return item.text
            elif isinstance(item, Cite):
                return ''.join(stringify(elem) for elem in item.content)
            elif str(item) == 'Space':
                return ' '
            else:
                return str(item)

        def caption():
            if elem.caption and elem.caption.content:
                return '\\caption{' + get_text(elem.caption.content) + '}\n' + \
                    '\\label{tbl:' + get_text(elem.caption.content).replace(' ', r'-').translate(str.maketrans("", "", "();,.")) + '}\n'
            return ''

        def tabular():
            return '\\begin{tabular}{' + 'l' * elem.cols + '}\n\\toprule\n'

        def headers():
            if elem.head and elem.head.content:
                return ' & '.join([get_text(cell.content) 
                                   for cell in elem.head.content[0].content]
                                  ) + '\\\\\n\\midrule\n'
            return ''

        def items():
            rows = []
            for body in elem.content:
                for row in body.content:
                    rows.append(' & '.join([get_text(cell.content) for cell in row.content]) + '\\\\')
            return '\n'.join(rows) + '\n'

        result = '\\begin{table*}[t]\n\\centering\n' + \
                 caption() + \
                 tabular() + \
                 headers() + \
                 items() + \
                 '\\bottomrule\n\\end{tabular}\n\\end{table*}'

        print("Table processed successfully", file=sys.stderr)
        return RawBlock(result, 'latex')
    except Exception as e:
        print(f"Error processing table: {str(e)}", file=sys.stderr)
        return elem


def prepare(doc):
    pass

def action(elem, doc):
    if doc.format != 'latex':
        return None

    if isinstance(elem, Table):
        print("Table found!", file=sys.stderr) 
        return replace_longtables_with_tabular(elem, doc)

    return None

def finalize(doc):
    pass

# keep this structure: see http://scorreia.com/software/panflute/guide.html
def main(doc=None):
    return run_filter(action,
                         prepare=prepare,
                         finalize=finalize,
                         doc=doc)


if __name__ == '__main__':
    main()
