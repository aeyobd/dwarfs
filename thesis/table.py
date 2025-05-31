import sys
import re

from panflute import *


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
    elif isinstance(item, Quoted):
        # Handle quoted text with proper quote characters
        quote = '"' if item.quote_type == 'DoubleQuote' else "'"
        content = ''.join(get_text(c) for c in item.content)
        return f'{quote}{content}{quote}'
    elif isinstance(item, RawInline):
        return item.text
    elif isinstance(item, Cite):
        return ''.join(stringify(elem) for elem in item.content)
    elif str(item) == 'Space':
        return ' '
    else:
        return str(item)



def table_caption(elem):
    if elem.caption and elem.caption.content:
        raw_caption = get_text(elem.caption.content)
    else:
        return ""

    # Extract components using regex
    label = None
    short = None
    caption_text = raw_caption
    labels = re.findall(r'\\label\{([^}]+)\}', raw_caption)

    # Extract short caption (e.g., short="Short caption")
    if elem.caption.short_caption:
        short = stringify(elem.caption.short_caption)

    if len(labels) > 1:
        debug("multiple labels")

    if len(labels):
        label = labels[0]
        # should only replace once if only 1 label
        caption_text = caption_text.replace(f"\\label{{{label}}}", "")
    elif elem.identifier:
        debug("using ID instead of label")
        label = get_text(elem.identifier)
    else:
        debug("no label found")

    # Build LaTeX components
    parts = []
    if short:
        parts.append(f'\\caption[{short}]{{{caption_text}}}')
    else:
        parts.append(f'\\caption{{{caption_text}}}')
    
    if label:
        parts.append(f'\\label{{{label}}}')

    return '\n'.join(parts) + '\n'
   
    return '\n'.join(parts) + '\n'
    return ''



def table_headers(elem):
    if elem.head and elem.head.content:
        return ' & '.join([get_text(cell.content) 
                           for cell in elem.head.content[0].content]
                          ) + '\\\\\n\\midrule\n'
    return ''


def table_items(elem):
    rows = []
    for body in elem.content:
        for row in body.content:
            rows.append(' & '.join([get_text(cell.content) for cell in row.content]) + '\\\\')
    return '\n'.join(rows) + '\n'



def tabular_begin(elem):
    return '\\begin{tabular}{' + 'l' * elem.cols + '}\n\\toprule\n'


def replace_longtables_with_tabular(elem, doc):
    try:
        result = '\\begin{table*}[t]\n\\centering\n' + \
                 table_caption(elem) + \
                 tabular_begin(elem) + \
                 table_headers(elem) + \
                 table_items(elem) + \
                 '\\bottomrule\n\\end{tabular}\n\\end{table*}'

        print("Table processed successfully", file=sys.stderr)
        return RawBlock(result, 'latex')

    except Exception as e:
        print(f"Error processing table: {str(e)}", file=sys.stderr)
        return elem



def action(elem, doc):
    if doc.format != 'latex':
        return None

    if isinstance(elem, Table):
        print("Table found!", file=sys.stderr) 
        return replace_longtables_with_tabular(elem, doc)

    return None



def main(doc=None):
    return run_filter(action,
                         doc=doc)


if __name__ == '__main__':
    main()
