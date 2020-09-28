cwlVersion: v1.1
class: CommandLineTool
label: Apply base quality score recalibration to a BAM file
doc: |-
  Apply base quality score recalibration to a BAM file using GATK 3.8-1.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: apply_bqsr_bam
    dockerFile: |-
      FROM broadinstitute/gatk3:3.8-1
      
      RUN apt-get install samtools
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  fasta_file:
    type: File
    doc: |-
      FASTA file for reference genome.
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
  bqsr_table_file:
    type: File
    doc: |-
      File with BQSR table.
    format: edam:format_1964
  bam_file:
    type: File
    doc: |-
      The BAM file to be analyzed.
    format: edam:format_2572
    secondaryFiles:
      - .bai
  threads:
    type: int
    doc: |-
      The number of threads that GATK should use.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created.
      #Output_File=edam:format_2572
arguments:
  - shellQuote: false
    valueFrom: |-
      # There must be a better way around this, but this tool looks for the
      
      # sequence dictionary file without .fa in the name. This is a workaround.
      
      ln -s "$(inputs.fasta_file.path)" .
      
      ln -s "$(inputs.fasta_file.path)".fai .
      
      # This line uses a JavaScript expression.
      
      ln -s "$(inputs.fasta_file.path).dict" "$(inputs.fasta_file.basename.replace('.fa', '.dict'))"
      
      java -Xms128m -Xmx2g -jar /usr/GenomeAnalysisTK.jar -T PrintReads -R "$(inputs.fasta_file.basename)" -BQSR "$(inputs.bqsr_table_file.path)" -I "$(inputs.bam_file.path)" -o "$(inputs.output_file_name)" -nct $(inputs.threads)
      
      # We use samtools rather than sambamba because we can install it more easily in this image.
      
      samtools index "$(inputs.output_file_name)"
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
    format: edam:format_2572
  regular_output_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name).bai"
    doc: |-
      An index file.
    format: edam:format_3327
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: apply_bqsr_bam_output.txt
stderr: apply_bqsr_bam_error.txt
 
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
