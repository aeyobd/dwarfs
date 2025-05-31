import panflute as pf
import re


def extract_attributes(inlines):
    # Reverse search for attribute block
    attr_text = []
    indices = []
    
    # Find closing '}'
    end_pos = None
    for i in reversed(range(len(inlines))):
        if isinstance(inlines[i], pf.Str) and '}' in inlines[i].text:
            end_pos = i
            break
    
    if end_pos is None:
        return None, None, inlines
    
    # Find opening '{'
    start_pos = None
    for i in reversed(range(end_pos + 1)):
        if isinstance(inlines[i], pf.Str) and '{' in inlines[i].text:
            start_pos = i
            break

    if start_pos is None:
        return None, None, inlines
    

    not_id_pos = start_pos
    for i in range(start_pos, end_pos):
        if isinstance(inlines[i], pf.Space):
            not_id_pos = i
            break

    # Extract attribute content
    attr_elements = inlines[not_id_pos:end_pos]
    id_str = pf.stringify(inlines[start_pos:not_id_pos])[2:]
    remaining_inlines = pf.ListContainer(*inlines[:start_pos],
                                      *inlines[end_pos+1:])
    
    attrs_dict = parse_attributes(attr_elements)

    return id_str, attrs_dict, remaining_inlines


def parse_attributes(attr_elements):
    if len(attr_elements) == 0:
        return {}
    # Convert attribute elements to text
    splits = [0]

    for i in range(len(attr_elements)):
        if isinstance(attr_elements[i], pf.Space):
            splits.append(i)
    splits.append(len(attr_elements))

    attrs_dict = {}
    for i in range(len(splits)-1):
        attr_str = pf.stringify(attr_elements[splits[i]:splits[i+1]])
        if attr_str:
            key, val = attr_str.split("=")
            if val[0] == '"' and val[-1] == '"':
                val = val[1:-1]
            elif val[0] == "'" and val[-1] == "'":
                val = val[1:-1]
            else:
                val = val.strip()

            attrs_dict[key.strip()] = val

    return attrs_dict




def process_table(element, doc):
    if isinstance(element, pf.Table) and element.caption:
        caption = element.caption
        if caption.content:
            #pf.debug("table caption: ")
            # Get first paragraph of caption
            if isinstance(caption.content[0], pf.Plain):
                inlines = caption.content[0].content
            
                #pf.debug(caption)
                # Extract attributes from inlines
                id_str, attrs, new_inlines = extract_attributes(inlines)
                #pf.debug("attributes", attrs)
                #pf.debug("new inlines", pf.stringify(new_inlines))

                if attrs:
                    # Set short caption
                    if 'short' in attrs.keys():
                        short_text = attrs['short']
                        #pf.debug("short caption", short_text)
                    
                        element.caption = pf.Caption(pf.Plain(*new_inlines),
                                   short_caption=[pf.Str(short_text)],
                                                     )
                        element.identifier = id_str
                        #pf.debug(id_str)
                        #pf.debug("final : ", (element.caption))
                        #pf.debug("id : ", (element.identifier))
                #pf.debug()

    return element


def main(doc=None):
    r = pf.run_filter(process_table, doc=doc)
    return r


if __name__ == "__main__":
    main()
