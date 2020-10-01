cwlVersion: v1.2
class: CommandLineTool
label: Sort and index a BAM file
doc: |-
  Sort and index a BAM file.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: sort_bam
    dockerFile: |-
      FROM quay.io/biocontainers/sambamba:0.7.1--h148d290_2
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  bam_file:
    type: File
    doc: |-
      BAM file to be sorted.
    format: edam:format_2572
  threads:
    type: int
    doc: |-
      The number of threads that samtools should use when sorting.
  output_file_name:
    type: string
    doc: |-
      Name of the sorted BAM file that will be created.
      #Output_File=edam:format_2572
arguments:
  - shellQuote: false
    valueFrom: |-
      sambamba sort -t $(inputs.threads) -o "$(inputs.output_file_name)" "$(inputs.bam_file.path)"
      
      sambamba index -t $(inputs.threads) "$(inputs.output_file_name)"
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
    format: edam:format_2572
  regular_output_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name).bai"
    doc: |-
      An index file.
    format: edam:format_1964
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: sort_bam_output.txt
stderr: sort_bam_error.txt
 
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
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
