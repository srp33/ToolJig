cwlVersion: v1.2
class: CommandLineTool
label: Apply GATK BaseRecalibrator to a BAM file
doc: |-
  Detect systematic errors in base quality scores using the Genome Analysis Toolkit (GATK) BaseRecalibrator tool. We use GATK3 because it supports multi-threading in a way that's simpler to use.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: calculate_bqsr_table
    dockerFile: |-
      FROM broadinstitute/gatk3:3.8-1
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  bam_file:
    type: File
    format: edam:format_2572
    secondaryFiles:
      - .bai
    doc: |-
      The BAM file to be analyzed.
  fasta_file:
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
    doc: |-
      FASTA file for reference genome.
  known_sites_vcf_file_1:
    type: File
    format: edam:format_3016
    secondaryFiles:
      - .idx
    doc: |-
      A file that contains known variants sites. It will be used for recalibrating BAM files.
  known_sites_vcf_file_2:
    type: File
    format: edam:format_3016
    secondaryFiles:
      - .idx
    doc: |-
      A file that contains known variants sites. It will be used for recalibrating BAM files.
  known_sites_vcf_file_3:
    type: File
    format: edam:format_3016
    secondaryFiles:
      - .idx
    doc: |-
      A file that contains known variants sites. It will be used for recalibrating BAM files.
  output_file:
    type: string
    doc: |-
      Name of the output file that will be created.
      #Output_File=format: format: edam:format_1964
  threads:
    type: int
    doc: |-
      The number of threads that GATK should use.
arguments:
  - shellQuote: false
    valueFrom: |-
      for KS_FILE in $(inputs.known_sites_vcf_file_1.path) $(inputs.known_sites_vcf_file_2.path) $(inputs.known_sites_vcf_file_3.path); do KNOWN_SITES_ARG="$KNOWN_SITES_ARG -knownSites $KS_FILE"; done
      
      # There must be a better way around this, but BaseRecalibrator looks for the
      
      # sequence dictionary file without .fa in the name. This is a workaround.
      
      ln -s "$(inputs.fasta_file.path)" .
      
      ln -s "$(inputs.fasta_file.path)".fai .
      
      # This line uses a JavaScript expression.
      
      ln -s "$(inputs.fasta_file.path).dict" "$(inputs.fasta_file.basename.replace('.fa', '.dict'))"
      
      java -Xms128m -Xmx8g -jar /usr/GenomeAnalysisTK.jar -T BaseRecalibrator -R "$(inputs.fasta_file.basename)" -I "$(inputs.bam_file.path)" -o "$(inputs.output_file)" -nct $(inputs.threads) $KNOWN_SITES_ARG
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file)"
    doc: |-
      Output file matching the name specified in the "output_file" input.
    format: format: edam:format_1964
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: calculate_bqsr_bam_output.txt
stderr: calculate_bqsr_bam_error.txt
 
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
 
s:dateCreated: "2021-04-12"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl