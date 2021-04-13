cwlVersion: v1.2
class: Workflow
id: secondary_sqrt_workflow
label: Square root considers secondary files
doc: |-
  Calculates the square root of a number stored in a file and saves the result to an output file. It does the same for two secondary files. It then sums those values and writes the sum to a file. This demonstrates using secondary files within a workflow.
inputs:
  - id: number_file
    type: File
    format: edam:format_1964
    secondaryFiles:
      - .a
      - .b
  - id: output_file
    type: string
outputs:
  - id: sum_files__output_file
    type: File
    outputSource: sum_files/output_file
steps:
  calc_sqrt:
    run: secondary_sqrt_tool.cwl
    in:
      - id: number_file
        source: number_file
      - id: output_file
        default: random__kkRlosRALc
    out:
      [output_file]
  sum_files:
    run: sum_files_tool.cwl
    in:
      - id: number_file
        source: calc_sqrt/output_file
      - id: output_file
        source: output_file
    out:
      [output_file]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
s:dateCreated: "2021-03-29"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
