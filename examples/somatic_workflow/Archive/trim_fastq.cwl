cwlVersion: v1.2
class: CommandLineTool
label: Trim paired FASTQ files
doc: |-
  Download the FASTQ files from an online repository and trim adapter sequences and low-quality bases from FASTQ files using the atropos software. Only paired-end reads are supported. atropos: https://github.com/jdidion/atropos
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: trim_fastq
    dockerFile: |-
      FROM quay.io/biocontainers/atropos:1.1.25--py36h516909a_0
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  fastq_file_1:
    type: File
    doc: |-
      FASTQ file for the first ends of paired-end reads.
    format: edam:format_1931
  fastq_file_2:
    type: File
    doc: |-
      FASTQ file for the second ends of paired-end reads.
    format: edam:format_1931
  args:
    type: string
    doc: |-
      Miscellaneous arguments to be forwarded to the atropos software. Consult the atropos documentation for information about arguments that can be specified. Example: "-q 20 -Q 20 --minimum-length 40 --threads 4".
arguments:
  - shellQuote: false
    valueFrom: |-
      atropos $(inputs.args) -pe1 $(inputs.fastq_file_1.path) -pe2 $(inputs.fastq_file_2.path) -o $(inputs.fastq_file_1.basename) -p $(inputs.fastq_file_2.basename)
outputs:
  fastq_file_1:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file_1.basename)"
    doc: |-
      Trimmed FASTQ file for the first ends of paired-end reads.
    format: edam:format_1931
  fastq_file_2:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file_2.basename)"
    doc: |-
      Trimmed FASTQ file for the second ends of paired-end reads.
    format: edam:format_1931
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: trim_fastq_output.txt
stderr: trim_fastq_error.txt
 
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
