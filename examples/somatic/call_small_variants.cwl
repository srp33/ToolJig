cwlVersion: v1.2
class: CommandLineTool
label: Calls somatic variants using Mutect2 (GATK) for a tumor/normal pair
doc: |-
  Call somatic variants using Mutect2 (GATK) for a tumor/normal pair.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
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
    doc: |-
      FASTA file containing reference genome.
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
  normal_bam_file:
    type: File
    doc: |-
      BAM file for normal sample.
    format: edam:format_2572
    secondaryFiles:
      - .bai
  tumor_bam_file:
    type: File
    doc: |-
      BAM file for tumor sample.
    format: edam:format_2572
    secondaryFiles:
      - .bai
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
  output_file:
    type: string
    doc: |-
      Name of the output file that will be created.
      #Output_File=edam:format_3016
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
      
      gatk FilterMutectCalls -R "$(inputs.fasta_file.basename)" -V unfiltered.vcf -O "$(inputs.output_file)"
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file)"
    doc: |-
      Output file matching the name specified in the "output_file" input.
    format: edam:format_3016
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: call_small_variants_output.txt
stderr: call_small_variants_error.txt
 
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
 
s:dateCreated: "2020-07-13"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
