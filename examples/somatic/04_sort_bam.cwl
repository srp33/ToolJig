cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Sort a BAM file.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 04_sort_bam
    dockerFile: |-
      FROM quay.io/biocontainers/samtools:1.3--h0592bc0_3
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  bam_file:
    type: File
    doc: |-
      BAM file to be sorted.
  threads:
    type: int
    doc: |-
      The number of threads that samtools should use when sorting.
  output_file_name:
    type: string
    doc: |-
      Name of the BAM file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        samtools sort -@ $(inputs.threads) -o "$(inputs.output_file_name)" "$(inputs.bam_file.path)"
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Here we indicate that an output file matching the name specified in the inputs should be generated.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: 04_output.txt
stderr: 04_error.txt
