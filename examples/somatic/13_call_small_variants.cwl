cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Calls somatic variants using Mutect2 (GATK) for a tumor/normal pair.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 13_call_somatic_variants
    dockerFile: |-
      FROM broadinstitute/gatk:4.1.6.0
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  ref_genome_dir:
    type: Directory
    doc: |-
      Directory containing reference genome FASTA file and index files.
  ref_genome_fasta_name:
    type: string
    doc: |-
      Name of the FASTA file containing the reference genome.
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
        gatk Mutect2 -R "$(inputs.ref_genome_dir.path)/$(inputs.ref_genome_fasta_name)" --input $(inputs.normal_bam_file.path) --input $(inputs.tumor_bam_file.path) --output unfiltered.vcf -normal $(inputs.normal_sample_id) -tumor $(inputs.tumor_sample_id) --native-pair-hmm-threads $(inputs.threads)

        gatk FilterMutectCalls -R "$(inputs.ref_genome_dir.path)/$(inputs.ref_genome_fasta_name)" -V unfiltered.vcf -O "$(inputs.output_file_name)"
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
