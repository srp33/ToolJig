cwlVersion: v1.1
class: Workflow
id: simple_workflow_example
label: Simple workflow example
doc: |-
  This workflow accepts four integers, adds the first two, adds the last two, and then creates a text file that stores the two sums.
inputs:
  add1__number1: int
  add1__number2: int
  add1__output_file_name: string
  add2__number1: int
  add2__number2: int
  add2__output_file_name: string
  cat__output_file_name: string
outputs:
  cat__output_from_input_1:
    type: File
    outputSource: cat/output_from_input_1
steps:
  add1:
    run: add_tool.cwl
    in:
      number1: add1__number1
      number2: add1__number2
      output_file_name: add1__output_file_name
    out:
      [output_from_input_1]
  add2:
    run: add_tool.cwl
    in:
      number1: add2__number1
      number2: add2__number2
      output_file_name: add2__output_file_name
    out:
      [output_from_input_1]
  cat:
    run: cat_tool.cwl
    in:
      file1: add1/output_from_input_1
      file2: add2/output_from_input_1
      output_file_name: cat__output_file_name
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