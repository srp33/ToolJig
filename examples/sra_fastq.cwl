cwlVersion: v1.1
class: CommandLineTool
doc: |-
  This tool downloads an archive from the Sequence Read Archive and converts it to a set of FASTQ files.
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
arguments:
  - shellQuote: false
    valueFrom: |-
      fasterq-dump --split-files $(inputs.srr_id)

      gzip *.fastq
outputs:
  output_1:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.fastq.gz"
    doc: |-
      Output files matching the value specified in the "output_file_pattern" input should be generated.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: output.txt
stderr: error.txt
