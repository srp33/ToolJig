#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: Workflow
id: simple_workflow_example
label: Simple workflow example
doc: |-
  This workflow accepts four integers, adds the first two, adds the last two, and then creates a text file that contains the two sums.
inputs:
  number1: int
  number2: int
  number3: int
  number4: int
  sum1_file_name: string
  sum2_file_name: string
  output_file_name: string
outputs:
  cat_result:
    type: File
    outputSource: cat/output_from_input_1
steps:
  add1:
    run: add_tool.cwl
    in:
      number1: number1
      number2: number2
      output_file_name: sum1_file_name
    out: [output_from_input_1]
  add2:
    run: add_tool.cwl
    in:
      number1: number3
      number2: number4
      output_file_name: sum2_file_name
    out: [output_from_input_1]
  cat:
    run: cat_tool.cwl
    in:
      file1: add1/output_from_input_1
      file2: add2/output_from_input_1
      output_file_name: output_file_name
    out: [output_from_input_1]