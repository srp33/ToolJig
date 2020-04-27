cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Mark duplicate reads in a BAM file. It also indexes the file.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 06_mark_dups
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
  threads:
    type: int
    doc: |-
      The number of threads that sambamba will use.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        sambamba markdup -t $(inputs.threads) $(inputs.bam_file.path) "$(inputs.output_file_name)"

        sambamba index -t $(inputs.threads) "$(inputs.output_file_name)" 
outputs:
  output_file_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Duplicate-marked BAM file.
  output_file_2:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name).bai"
    doc: |-
      Index of duplicate-marked BAM file.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: output.txt
stderr: error.txt
