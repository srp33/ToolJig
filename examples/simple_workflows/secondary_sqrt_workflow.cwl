cwlVersion: v1.2
class: Workflow
id: secondary_sqrt_workflow
label: Square root considers secondary files
doc: |-
  Calculates the square root of a number stored in a file and saves the result to an output file. It does the same for two secondary files. This demonstrates using secondary files from a tool within a workflow.
inputs:
  - id: calculation1__number_file
    type: File
    format: edam:format_1964
    secondaryFiles:
      - .a
      - .b
  - id: calculation1__output_file
    type: string
outputs:
  - id: calculation1__output_file
    type: File
    outputSource: calculation1/output_file
steps:
  calculation1:
    run: secondary_sqrt_tool.cwl
    in:
      - id: number_file
        source: calculation1__number_file
      - id: output_file
        source: calculation1__output_file
    out:
      [output_file]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
s:dateCreated: "2021-02-04"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl