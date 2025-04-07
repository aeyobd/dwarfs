from panflute import *
import sys
from os import path

def rename_url(blocks, i):
    content = blocks[i].content
    for j in range(len(content)):
        if isinstance(content[j], Plain):
            for con2 in content[j].content:
                if isinstance(con2, Image):
                    url = con2.url
                    print("Found image at ", url, file=sys.stderr)
                    if url.startswith("figures"):
                        url_new = path.splitext(url)[0] + ".pdf"
                        if path.isfile(url_new):
                            con2.url = url_new
                            print("Renaming ", url,"=>", con2.url, file=sys.stderr)


def process_figures(doc):
    blocks = list(doc.content)
    i = 0
    while i < len(blocks):
        if isinstance(blocks[i], Figure):
            print("Figure found!", file=sys.stderr) 
            url = rename_url(blocks, i)

            if i + 1 < len(blocks) and isinstance(blocks[i+1], Para):
                next_para = blocks[i+1]
                if next_para.content:
                    first_part = next_para.content[0]
                    if (len(next_para.content) > 2 
                            and isinstance(first_part, Str) 
                            and first_part.text == "Figure:"
                            and isinstance(next_para.content[1], Space)
                        ):

                        print("Figure caption found!", file=sys.stderr) 
                         # Get original caption from the Figure
                        original_caption = blocks[i].caption
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
