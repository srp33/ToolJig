cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Download a file from an online location.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: download_fastq
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  url:
    type: string
    doc: |-
      URL of file to be download.
  out_file_name:
    type: string
    doc: |-
      Name of the output file to be saved.
arguments:
  - shellQuote: false
    valueFrom: |-
      wget -O $(inputs.out_file_name) $(inputs.url)
outputs:
  output_1:
    type: File
    outputBinding:
      glob: "$(inputs.out_file_name)"
    doc: |-
      Output file matching the name specified in the "out_file_name" input.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: output.txt
stderr: error.txt
