cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Prepare a reference genome for variant calling.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 02_prep_ref_genome
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
  output_file:
    type: 
      type: array
      items: File
    outputBinding:
      glob: "$(inputs.ref_genome_version).fa*"
    doc: |-
      FASTA file
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: 02_output.txt
stderr: 02_error.txt
