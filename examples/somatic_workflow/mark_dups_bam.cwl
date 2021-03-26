cwlVersion: v1.2
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
    format: edam:format_2572
    doc: |-
      The BAM file to be marked.
  output_file:
    type: string
    doc: |-
      Name of the output file that will be created.
      #Output_File=format: edam:format_2572;secondaryFiles: .bai
  threads:
    type: int
    doc: |-
      The number of threads that sambamba will use.
arguments:
  - shellQuote: false
    valueFrom: |-
      sambamba markdup -t $(inputs.threads) $(inputs.bam_file.path) "$(inputs.output_file)"
      
      sambamba index -t $(inputs.threads) "$(inputs.output_file)" 
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file)"
    doc: |-
      Output file matching the name specified in the "output_file" input.
    format: edam:format_2572
    secondaryFiles:
      - .bai
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
 
s:dateCreated: "2021-03-03"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl