cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Index a BAM file.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 05_index_bam
    dockerFile: |-
      FROM quay.io/biocontainers/sambamba:0.7.1--h148d290_2
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  bam_file:
    type: File
    doc: |-
      BAM file to be indexed.
  threads:
    type: int
    doc: |-
      The number of threads to be used by samtools.
  output_file_name:
    type: string
    doc: |-
      Name of the index file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        sambamba index -t $(inputs.threads) "$(inputs.bam_file.path)" "$(inputs.output_file_name)"
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
stdout: output.txt
stderr: error.txt
