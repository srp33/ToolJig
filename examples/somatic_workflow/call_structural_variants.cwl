cwlVersion: v1.2
class: CommandLineTool
label: Call structural variants using Delly for a tumor/normal pair
doc: |-
  Call structural variants using Delly for a tumor/normal pair.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: call_structural_variants
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
      
      RUN conda install -c bioconda/label/cf201901 delly bcftools -y
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
  exclude_template_url:
    type: string
    doc: |-
      URL of a file with regions to exclude (which speeds up execution times). These files can be found for many reference genomes here: https://github.com/dellytools/delly/tree/master/excludeTemplates. Make sure to get the URL for the raw version of the file. Example: https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created.
      #Output_File=edam:format_3016
arguments:
  - shellQuote: false
    valueFrom: |-
      wget $(inputs.exclude_template_url)
      
      delly call -x `basename "$(inputs.exclude_template_url)"` -o output.bcf -g "$(inputs.fasta_file.path)" "$(inputs.tumor_bam_file.path)" "$(inputs.normal_bam_file.path)"
      
      bcftools view output.bcf > "$(inputs.output_file_name)"
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
    format: edam:format_3016
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: call_structural_variants_output.txt
stderr: call_structural_variants_error.txt
 
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
