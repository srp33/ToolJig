cwlVersion: v1.2
class: Workflow
id: recursive_sqrt_root_workflow
label: Recursive square root
doc: |-
  This workflow reads a number from a file, calculates the square root of that number, calculates the square root of the resulting number, and saves the output to a file. This demonstrates the ability to invoke the same tool recursively.
inputs:
  - id: number_file
    type: File
    format: edam:format_1964
  - id: output_file
    type: string
outputs:
  - id: calculation2__output_file
    type: File
    outputSource: calculation2/output_file
steps:
  calculation1:
    run: sqrt_tool.cwl
    in:
      - id: number_file
        source: number_file
      - id: output_file
        default: random__RmAxPgNnvI
    out:
      [output_file]
  calculation2:
    run: sqrt_tool.cwl
    in:
      - id: number_file
        source: calculation1/output_file
      - id: output_file
        source: output_file
    out:
      [output_file]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
s:dateCreated: "2021-03-30"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
