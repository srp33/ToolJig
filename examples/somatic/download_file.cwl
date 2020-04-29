cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Download a file from an online location.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: download_file
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
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: download_file_output.txt
stderr: download_file_error.txt
