cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Calls somatic variants using Mutect2 (GATK) for a tumor/normal pair.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: call_small_variants
    dockerFile: |-
      FROM broadinstitute/gatk:4.1.6.0
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  fasta_file:
    type: File
    secondaryFiles:
      - .fai
      - .dict
    doc: |-
      FASTA file containing reference genome.
  normal_bam_file:
    type: File
    secondaryFiles:
      - .bai
    doc: |-
      BAM file for normal sample.
  tumor_bam_file:
    type: File
    secondaryFiles:
      - .bai
    doc: |-
      BAM file for tumor sample.
  normal_sample_id:
    type: string
    doc: |-
      Unique identifier for the normal sample.
  tumor_sample_id:
    type: string
    doc: |-
      Unique identifier for the tumor sample.
  threads:
    type: int
    doc: |-
      The number of threads that GATK should use.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        # There must be a better way around this, but this tool looks for the
        # sequence dictionary file without .fa in the name. This is a workaround.
        ln -s "$(inputs.fasta_file.path)" .

        ln -s "$(inputs.fasta_file.path)".fai .

        # This line uses a JavaScript expression.
        ln -s "$(inputs.fasta_file.path).dict" "$(inputs.fasta_file.basename.replace('.fa', '.dict'))"

        gatk Mutect2 -R "$(inputs.fasta_file.basename)" --input $(inputs.normal_bam_file.path) --input $(inputs.tumor_bam_file.path) --output unfiltered.vcf -normal $(inputs.normal_sample_id) -tumor $(inputs.tumor_sample_id) --native-pair-hmm-threads $(inputs.threads)

        gatk FilterMutectCalls -R "$(inputs.fasta_file.basename)" -V unfiltered.vcf -O "$(inputs.output_file_name)"
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: call_small_variants_output.txt
stderr: call_small_variants_error.txt
