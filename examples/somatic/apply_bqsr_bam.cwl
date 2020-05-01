cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Apply base quality score recalibration to a BAM file.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
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
    secondaryFiles:
      - .fai
      - .dict
    doc: |-
      FASTA file for reference genome.
  bqsr_table_file:
    type: File
    doc: |-
      File with BQSR table.
  bam_file:
    type: File
    secondaryFiles:
      - .bai
    doc: |-
      The BAM file to be analyzed.
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

        java -Xms128m -Xmx2g -jar /usr/GenomeAnalysisTK.jar -T PrintReads -R "$(inputs.fasta_file.basename)" -BQSR "$(inputs.bqsr_table_file.path)" -I "$(inputs.bam_file.path)" -o "$(inputs.output_file_name)" -nct $(inputs.threads)

        # We use samtools rather than sambamba because we can install it more easily in this image.
        samtools index "$(inputs.output_file_name)"
outputs:
  output_file_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    secondaryFiles:
      - .bai
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: apply_bqsr_bam_output.txt
stderr: apply_bqsr_bam_error.txt
