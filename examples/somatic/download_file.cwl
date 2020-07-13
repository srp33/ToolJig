cwlVersion: v1.1
class: CommandLineTool
label: Download a file from an online location
doc: |-
  Download a file from an online location.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: download_file
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  url:
    type: string
    doc: |-
      URL of file to be download.
  out_file_name:
    type: string
    doc: |-
      Name of the output file to be saved.
      #Output_File=edam:format_1964
arguments:
  - shellQuote: false
    valueFrom: |-
      wget -O $(inputs.out_file_name) $(inputs.url)
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.out_file_name)"
    doc: |-
      Output file matching the name specified in the "out_file_name" input.
    format: edam:format_1964
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: download_file_output.txt
stderr: download_file_error.txt
 
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
 
s:dateCreated: "2020-07-13"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.23.owl