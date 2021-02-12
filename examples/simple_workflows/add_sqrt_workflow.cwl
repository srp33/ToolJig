cwlVersion: v1.2
class: Workflow
id: add_sqrt
label: Add two numbers, then take square root
doc: |-
  This workflow accepts two integers, adds them, calculates the square root of that number, and then stores the square root in a file.
inputs:
  - id: add__number1
    type: int
  - id: add__number2
    type: int
  - id: add__output_file
    type: string
  - id: sqrt__output_file
    type: string
outputs:
  - id: sqrt__output_file
    type: File
    outputSource: sqrt/output_file
steps:
  add:
    run: add_tool.cwl
    in:
      number1: add__number1
      number2: add__number2
      output_file: add__output_file
    out:
      [output_file]
  sqrt:
    run: sqrt_tool.cwl
    in:
      number_file: add/output_file
      output_file: sqrt__output_file
    out:
      [output_file]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
s:dateCreated: "2021-02-04"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl