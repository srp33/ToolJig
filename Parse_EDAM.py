import requests

url = "http://edamontology.org/EDAM_1.23.tsv"

edam_content = requests.get(url, allow_redirects=True).content.decode()

header_items = edam_content.split("\n")[0].split("\t")
obsolete_index = header_items.index("Obsolete")
prefix_index = header_items.index("http://data.bioontology.org/metadata/prefixIRI")
label_index = header_items.index("Preferred Label")
def_index = header_items.index("Definitions")

out_dict = {}

for line in edam_content.split("\n")[1:]:
    line_items = line.split("\t")
    if len(line_items) < 2:
        continue

    is_obsolete = line_items[obsolete_index] == "TRUE"
    the_format = line_items[prefix_index]
    is_format = the_format.startswith("format_")

    if not is_format or is_obsolete:
        continue

    label = line_items[label_index]
    definition = line_items[def_index].replace("\"", "").split("|")[0]

    if len(definition) > 77:
        definition = definition[:77]
        definition = " ".join(definition.split(" ")[:-1]) + "..."

    description = f"{label} - {definition}"

    out_dict[label] = f"                        <option value=\"edam:{the_format}\">{description}</option>"

out_lines = [out_dict[label] for label in sorted(out_dict, key=str.casefold)]

out_lines.insert(0, out_dict.pop("plain text format (unformatted)"))
out_lines.insert(0, f"                        <option value=\"\"></option>")

for line in out_lines:
    print(line)
