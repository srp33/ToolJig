cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Prepare a reference genome for read alignment and variant calling. This tool also prepares index files and a dictionary for the reference sequence.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: prep_ref_genome
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
      RUN conda install -c bioconda/label/cf201901 bwa samtools picard -y
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  ref_genome_version:
    type: string
    doc: |-
      A version identifier for a reference genome from the UCSC repository. Options include hg19, hg38, etc. See https://hgdownload.soe.ucsc.edu/downloads.html.
arguments:
  - shellQuote: false
    valueFrom: |-
      wget "http://hgdownload.cse.ucsc.edu/goldenPath/$(inputs.ref_genome_version)/bigZips/$(inputs.ref_genome_version).fa.gz"

      gunzip $(inputs.ref_genome_version).fa.gz

      bwa index -a bwtsw $(inputs.ref_genome_version).fa

      samtools faidx $(inputs.ref_genome_version).fa

      picard -Xms128m -Xmx2g CreateSequenceDictionary REFERENCE=$(inputs.ref_genome_version).fa OUTPUT=$(inputs.ref_genome_version).fa.dict
outputs:
  output_file_1:
    type: File
    outputBinding:
      glob: "$(inputs.ref_genome_version).fa"
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .fai
      - .pac
      - .sa
      - .dict
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: prep_ref_genome_output.txt
stderr: prep_ref_genome_error.txt
