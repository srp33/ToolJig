cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Download the FASTQ files from an online repository and trim adapter sequences and low-quality bases from FASTQ files using the atropos software. Only paired-end reads are supported. atropos: https://github.com/jdidion/atropos
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
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
      URL of FASTQ file for the first ends of paired-end reads.
  fastq_file_2:
    type: File
    doc: |-
      URL of FASTQ file for the second ends of paired-end reads.
  args:
    type: string
    doc: |-
      Miscellaneous arguments to be forwarded to the atropos software. Consult the atropos documentation for information about arguments that can be specified. Example: "-q 20 -Q 20 --minimum-length 40 --threads 4".
arguments:
  - shellQuote: false
    valueFrom: |-
      atropos $(inputs.args) -pe1 $(inputs.fastq_file_1.path) -pe2 $(inputs.fastq_file_2.path) -o $(inputs.fastq_file_1.basename) -p $(inputs.fastq_file_2.basename)
outputs:
  output_1:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file_1.basename)"
  output_2:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file_2.basename)"
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: trim_fastq_output.txt
stderr: trim_fastq_error.txt
