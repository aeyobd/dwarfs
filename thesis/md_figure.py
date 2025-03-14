from panflute import *
import sys

def process_figures(doc):
    blocks = list(doc.content)
    i = 0
    while i < len(blocks):
        if isinstance(blocks[i], Figure):
            print("Figure found!", file=sys.stderr) 
            if i + 1 < len(blocks) and isinstance(blocks[i+1], Para):
                next_para = blocks[i+1]
                if next_para.content:
                    first_part = next_para.content[0]
                    if isinstance(first_part, Str) and first_part.text == "Figure:" and isinstance(next_para.content[1], Space):

                        print("Figure caption found!", file=sys.stderr) 
                         # Get original caption from the Figure
                        original_caption = blocks[i].caption
                        print(original_caption.content, file=sys.stderr)
                        print(original_caption.content[0], file=sys.stderr)
                        print(type(original_caption.content[0]), file=sys.stderr)
                        # Create new caption: [original] + LineBreak + [new content]
                        inl = original_caption.content[0].content

                        new_caption = Caption(
                                Plain(*next_para.content[2:]),
                                short_caption = [x for x in inl]
                                )
                        blocks[i].caption = new_caption
                        
                        # Remove the caption paragraph
                        del blocks[i+1]

        i += 1
    doc.content = blocks

def action(elem, doc):
    if doc.format != 'latex':
        return None

    return None


def prepare(doc):
    pass

def finalize(doc):
    pass

def main(doc=None):
    return run_filter(action,
        prepare=prepare, 
        finalize=process_figures,
        doc=doc, 
    )

if __name__ == "__main__":
    main()
