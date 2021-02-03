cwlVersion: v1.2
class: Workflow
id: secondary_sqrt
label: Square root considers secondary files
doc: |-
  Calculates the square root of a number stored in a file plus numbers stored in two secondary files. This demonstrates using secondary files from a tool within a workflow.
inputs:
  - id: calculation1__number_file
    type: File
    format: edam:format_1964
    secondaryFiles:
      - .a
      - .b
  - id: calculation1__output_file_name
    type: string
outputs:
  calculation1__output_from_input_1:
    type: File
    outputSource: calculation1/output_from_input_1
steps:
  calculation1:
    run: secondary_sqrt_tool.cwl
    in:
      number_file: calculation1__number_file
      output_file_name: calculation1__output_file_name
    out:
      [output_from_input_1]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
s:dateCreated: "2021-01-25"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
