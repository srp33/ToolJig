cwlVersion: v1.2
class: CommandLineTool
label: Prepare a reference genome for read alignment and variant calling
doc: |-
  Prepare a reference genome for read alignment and variant calling. This tool also prepares index files and a dictionary for the reference sequence.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
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
  regular_output_1:
    type: File
    outputBinding:
      glob: "$(inputs.ref_genome_version).fa"
    doc: |-
      FASTA file with accompanying index files and dictionary.
    format: edam:format_1929
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
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: prep_ref_genome_output.txt
stderr: prep_ref_genome_error.txt
 
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
