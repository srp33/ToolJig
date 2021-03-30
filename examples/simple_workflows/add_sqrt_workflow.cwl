cwlVersion: v1.2
class: Workflow
id: add_sqrt_workflow
label: Add two numbers, then take square root
doc: |-
  This workflow accepts two integers, adds them, calculates the square root of that number, and then stores the square root in a file.
inputs:
  - id: number1
    type: int
  - id: number2
    type: int
  - id: output_file
    type: string
outputs:
  - id: sqrt__output_file
    type: File
    outputSource: sqrt/output_file
steps:
  add:
    run: add_tool.cwl
    in:
      - id: number1
        source: number1
      - id: number2
        source: number2
      - id: output_file
        default: temp_ZzutRAtezp
    out:
      [output_file]
  sqrt:
    run: sqrt_tool.cwl
    in:
      - id: number_file
        source: add/output_file
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
