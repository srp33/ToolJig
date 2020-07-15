cwlVersion: v1.1
class: CommandLineTool
label: Mark duplicate reads in a BAM file. It also indexes the file
doc: |-
  Mark duplicate reads in a BAM file. It also indexes the file.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: mark_dups_bam
    dockerFile: |-
      FROM quay.io/biocontainers/sambamba:0.7.1--h148d290_2
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  bam_file:
    type: File
    doc: |-
      The BAM file to be marked.
    format: edam:format_2572
  threads:
    type: int
    doc: |-
      The number of threads that sambamba will use.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created.
      #Output_File=edam:format_2572
arguments:
  - shellQuote: false
    valueFrom: |-
      sambamba markdup -t $(inputs.threads) $(inputs.bam_file.path) "$(inputs.output_file_name)"
      
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
    format: edam:format_3327
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: mark_dups_bam_output.txt
stderr: mark_dups_bam_error.txt
 
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
