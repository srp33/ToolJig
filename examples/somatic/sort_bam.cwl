cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Sort and index a BAM file.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
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
  threads:
    type: int
    doc: |-
      The number of threads that samtools should use when sorting.
  output_file_name:
    type: string
    doc: |-
      Name of the sorted BAM file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        sambamba sort -t $(inputs.threads) -o "$(inputs.output_file_name)" "$(inputs.bam_file.path)"

        sambamba index -t $(inputs.threads) "$(inputs.output_file_name)"
outputs:
  output_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
  output_2:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name).bai"
    doc: |-
      An index file.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: sort_bam_output.txt
stderr: sort_bam_error.txt
