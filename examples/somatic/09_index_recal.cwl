cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Index a VCF file using GATK's IndexFeatureFile tool.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 09_index_recal
    dockerFile: |-
      FROM broadinstitute/gatk:4.1.6.0
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  vcf_file:
    type: File
    doc: |-
      VCF file that needs to be indexed.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        gatk IndexFeatureFile --input $(inputs.vcf_file.path) --output $(inputs.output_file_name)
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
stdout: 09_output.txt
stderr: 09_error.txt
