cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Download the FASTQ files from an online repository and trim adapter sequences and low-quality bases from FASTQ files using the atropos software. Only paired-end reads are supported. atropos: https://github.com/jdidion/atropos
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 01_trim_fastq
    dockerFile: |-
      FROM quay.io/biocontainers/atropos:1.1.25--py36h516909a_0
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  fastq_url_1:
    type: File
    doc: |-
      URL of FASTQ file for the first ends of paired-end reads.
  fastq_url_2:
    type: File
    doc: |-
      URL of FASTQ file for the second ends of paired-end reads.
  args:
    type: string
    doc: |-
      Miscellaneous arguments to be forwarded to the atropos software. Consult the atropos documentation for information about arguments that can be specified. Example: "-q 20 -Q 20 --minimum-length 40 --threads 4".
  out_fastq_file_name_1:
    type: string
    doc: |-
      Name of the first output (FASTQ) file with trimmed reads.
  out_fastq_file_name_2:
    type: string
    doc: |-
      Name of the second output (FASTQ) file with trimmed reads.
arguments:
  - shellQuote: false
    valueFrom: |-
      

      atropos $(inputs.args) -pe1 $(inputs.fastq_file_1.path) -pe2 $(inputs.fastq_file_2.path) -o $(inputs.out_fastq_file_name_1) -p $(inputs.out_fastq_file_name_2)
outputs:
  output_1:
    type: File
    outputBinding:
      glob: "$(inputs.out_fastq_file_name_1)"
    doc: |-
      Output file matching the name specified in the "out_fastq_file_name_1" input.
  output_2:
    type: File
    outputBinding:
      glob: "$(inputs.out_fastq_file_name_2)"
    doc: |-
      Output file matching the name specified in the "out_fastq_file_name_2" input.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: 01_output.txt
stderr: 01_error.txt
