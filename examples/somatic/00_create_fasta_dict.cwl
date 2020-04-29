cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Create a FASTA sequence dictionary using Picard tools.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: create_fasta_dict
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
      RUN conda install -c bioconda/label/cf201901 picard -y
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  fasta_file:
    type: File
    doc: |-
      FASTA file.
  output_file_name:
    type: string
    doc: |-
      Name of the output file.
arguments:
  - shellQuote: false
    valueFrom: |-
      picard -Xms128m -Xmx2g CreateSequenceDictionary REFERENCE=$(inputs.fasta_file.path) OUTPUT=$(inputs.output_file_name)
outputs:
  output_file_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: create_fasta_dict_output.txt
stderr: create_fasta_dict_error.txt
