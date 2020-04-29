cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Mark duplicate reads in a BAM file. It also indexes the file.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
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
    secondaryFiles:
      - .bai
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: mark_dups_bam_output.txt
stderr: mark_dups_bam_error.txt
