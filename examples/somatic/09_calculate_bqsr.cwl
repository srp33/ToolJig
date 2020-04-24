cwlVersion: v1.1
class: CommandLineTool
doc: Detect systematic errors in base quality scores using the Genome Analysis Toolkit (GATK) BaseRecalibrator tool. We use GATK3 because it supports multi-threading in a way that's simpler to take advantage of.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: 09_calculate_bqsr
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
      Path to a directory that contains VCF files with known polymorphic sites. Each file with .vcf extension will be used for recalibration.
  bam_file:
    type: File
    doc: |-
      The BAM file to be analyzed.
  threads:
    type: int
    doc: |
      The number of threads that GATK should use.
  args:
    type: string
    doc: |
      Additional arguments that will be passed through to the GATK (BaseRecalibrator) tool. Example value: "-L chr20" (limits recalibration to chromosome 20).
  output_file_name:
    type: string
    doc: |
      Name of the output file that will be created.
arguments:
    - shellQuote: false
      valueFrom: >
        KNOWN_SITES_COMMANDS=()

        for KS_FILE in ${KNOWN_SITES_FILE[@]}; do
    if [[ ! -f /volumes/vcf_dir/"${KS_FILE}.idx" ]]; then
        gatk3 -Xms128m -Xmx2g IndexFeatureFile -F /volumes/vcf_dir/"${KS_FILE}"
    fi
    KNOWN_SITES_COMMANDS+=("--knownSites /volumes/vcf_dir/${KS_FILE}")
done

gatk3 -Xms128m -Xmx${MAX_MEMORY} 

      -T BaseRecalibrator 

      $ARGS 

      -R /volumes/ref_dir/"${FASTA_FILE}" 

      -I /volumes/input_bam_dir/"${INPUT_BAM_FILE}" 

      ${KNOWN_SITES_COMMANDS[@]} 

      -o /volumes/output_bam_dir/"${OUTPUT_BAM_FILE}" 

      -nct ${THREADS}
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: Here we indicate that an output file matching the name specified in the inputs should be generated.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: 08_output.txt
stderr: 08_error.txt
