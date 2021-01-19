cwlVersion: v1.2
class: Workflow
id: simple_workflow_example
label: Simple workflow example
doc: |-
  This workflow accepts two integers, adds them, calculates the square root of that number, and then stores the square root in a file.
inputs:
  add__number1: int
  add__number2: int
  add__output_file_name: string
  sqrt__output_file_name: string
outputs:
  sqrt__output_from_input_1:
    type: File
    outputSource: sqrt/output_from_input_1
steps:
  add:
    run: add_tool.cwl
    in:
      number1: add__number1
      number2: add__number2
      output_file_name: add__output_file_name
    out:
      [output_from_input_1]
  sqrt:
    run: sqrt_tool.cwl
    in:
      number_file: add/output_from_input_1
      output_file_name: sqrt__output_file_name
    out:
      [output_from_input_1]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
s:dateCreated: "2020-08-05"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl