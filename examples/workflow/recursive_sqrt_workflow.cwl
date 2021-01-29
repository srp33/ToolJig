cwlVersion: v1.2
class: Workflow
id: recursive_square_root
label: Recursive square root
doc: |-
  This workflow reads a number from a file, calculates the square root of that number, calculates the square root of the resulting number, and saves the output to a file. This demonstrates the ability to invoke the same tool recursively.
inputs:
    - id: calculation1__number_file
      type: File
      format: edam:format_1964
    - id: calculation1__output_file_name
      type: string
    - id: calculation2__output_file_name
      type: string
outputs:
  calculation2__output_from_input_1:
    type: File
    outputSource: calculation2/output_from_input_1
steps:
  calculation1:
    run: sqrt_tool.cwl
    in:
      number_file: calculation1__number_file
      output_file_name: calculation1__output_file_name
    out:
      [output_from_input_1]
  calculation2:
    run: sqrt_tool.cwl
    in:
      number_file: calculation1/output_from_input_1
      output_file_name: calculation2__output_file_name
    out:
      [output_from_input_1]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
s:dateCreated: "2021-01-29"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl