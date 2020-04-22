cwlVersion: v1.1
class: CommandLineTool
doc: |-
  This tool downloads an archive from the Sequence Read Archive and converts it to a set of FASTQ files. It's designed for Illumina sequencing data processed as paired-end reads.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: sra_fastq
    dockerFile: |-
      FROM quay.io/biocontainers/sra-tools:2.10.3--pl526haddd2b5_0
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  srr_id:
    type: string
    doc: |-
      This is a unique identifier for a sample run. It should have an "SRR" prefix (e.g., SRR11180057).
  fastq_file_name_1:
    type: string
    doc: |-
      This is the name of a FASTQ file that will be created for the first end of paired-end reads. This file will be gzipped. Example: SRR11180057_1.fastq.gz.
  fastq_file_name_2:
    type: string
    doc: |-
      This is the name of a FASTQ file that will be created for the second end of paired-end reads. Example: SRR11180057_2.fastq.gz.
arguments:
  - shellQuote: false
    valueFrom: |-
      fasterq-dump --split-files $(inputs.srr_id)

      gzip *.fastq
outputs:
  output_1:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file_name_1)"
    doc: |-
      Output file matching the name specified in the "fastq_file_name_1" input.
  output_2:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file_name_2)"
    doc: |-
      Output file matching the name specified in the "fastq_file_name_2" input.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: 01_output.txt
stderr: 01_error.txt
