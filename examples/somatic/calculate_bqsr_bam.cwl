cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Detect systematic errors in base quality scores using the Genome Analysis Toolkit (GATK) BaseRecalibrator tool. We use GATK3 because it supports multi-threading in a way that's simpler to take advantage of.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: calculate_bqsr_bam
    dockerFile: |-
      FROM broadinstitute/gatk3:3.8-1
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

        java -Xms128m -Xmx2g -jar /usr/GenomeAnalysisTK.jar -T BaseRecalibrator -R "$(inputs.ref_genome_dir.path)/$(inputs.ref_genome_fasta_name)" -I "$(inputs.bam_file.path)" -o "$(inputs.output_file_name)" -nct $(inputs.threads) $KNOWN_SITES_ARG
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
stdout: calculate_bqsr_bam_output.txt
stderr: calculate_bqsr_bam_error.txt
