cwlVersion: v1.2
class: CommandLineTool
label: Download a FASTQ file from an online location
doc: |-
  Download a FASTQ file from an online location.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: download_fastq_file
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  output_file:
    type: string
    doc: |-
      Name of the output file to be saved.
      #Output_File=format: edam:format_1931
  url:
    type: string
    doc: |-
      URL of the file to be downloaded. This file must be in FASTQ format and must be compressed using the gzip utility.
arguments:
  - shellQuote: false
    valueFrom: |-
      wget -O $(inputs.output_file) $(inputs.url)
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file)"
    doc: |-
      Output file matching the name specified in the "output_file" input.
    format: edam:format_1931
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
 
s:dateCreated: "2021-03-02"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl