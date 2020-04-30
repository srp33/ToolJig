cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Detect systematic errors in base quality scores using the Genome Analysis Toolkit (GATK) BaseRecalibrator tool. We use GATK3 because it supports multi-threading in a way that's simpler to use.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: calculate_bqsr_bam
    dockerFile: |-
      FROM broadinstitute/gatk3:3.8-1
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
  known_sites_dir:
    type: Directory
    doc: |-
      Path to a directory that contains VCF files with known polymorphic sites. Each file with a .vcf extension will be used for recalibration.
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
        for KS_FILE in $(inputs.known_sites_dir.path)/*.vcf; do KNOWN_SITES_ARG="$KNOWN_SITES_ARG -knownSites $KS_FILE"; done

        # There must be a better way around this, but BaseRecalibrator looks for the
        # sequence dictionary file without .fa in the name. This is a workaround.
        ln -s "$(inputs.fasta_file.path)" .

        ln -s "$(inputs.fasta_file.path)".fai .

        # This line uses a JavaScript expression.
        ln -s "$(inputs.fasta_file.path).dict" "$(inputs.fasta_file.basename.replace('.fa', '.dict'))"

        java -Xms128m -Xmx8g -jar /usr/GenomeAnalysisTK.jar -T BaseRecalibrator -R "$(inputs.fasta_file.basename)" -I "$(inputs.bam_file.path)" -o "$(inputs.output_file_name)" -nct $(inputs.threads) $KNOWN_SITES_ARG
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: calculate_bqsr_bam_output.txt
stderr: calculate_bqsr_bam_error.txt
