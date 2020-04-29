cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Call structural variants using Delly for a tumor/normal pair.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: call_structural_variants
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
      RUN conda install -c bioconda/label/cf201901 delly bcftools -y
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
  exclude_template_url:
    type: string
    doc: |-
      URL of a file with regions to exclude (which speeds up execution times). These files can be found for many reference genomes here: https://github.com/dellytools/delly/tree/master/excludeTemplates. Make sure to get the URL for the raw version of the file. Example: https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        wget $(inputs.exclude_template_url)

        delly call -x `basename "$(inputs.exclude_template_url)"` -o output.bcf -g "$(inputs.ref_genome_dir.path)/$(inputs.ref_genome_fasta_name)" "$(inputs.tumor_bam_file.path)" "$(inputs.normal_bam_file.path)"

        bcftools view output.bcf > "$(inputs.output_file_name)"
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
stdout: call_structural_variants_output.txt
stderr: call_structural_variants_error.txt
