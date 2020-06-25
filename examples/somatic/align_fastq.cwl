cwlVersion: v1.1
class: CommandLineTool
doc: |-
  Align FASTQ files to a reference genome using the Burrows-Wheeler Aligner software (bwa mem). This is designed for paired-end reads stored in two separate FASTQ files.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: align_fastq
    dockerFile: |-
      FROM biocontainers/biocontainers:v1.0.0_cv4
      
      RUN conda install -c bioconda/label/cf201901 bwa samtools -y
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
inputs:
  fasta_file:
    type: File
    doc: |-
      Reference genome FASTA file.
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .fai
      - .pac
      - .sa
  fastq_file_1:
    type: File
    doc: |-
      The first FASTQ file.
  fastq_file_2:
    type: File
    doc: |-
      The second FASTQ file.
  threads:
    type: int
    doc: |-
      The number of threads that BWA should use during alignment.
  read_group_string:
    type: string
    doc: |-
      This argument allows you to specify read-group information during alignment. Please specify the whole read-group string, including the @RG prefix, surround in quotes. You can find a helpful tutorial here: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472.
  args:
    type: string
    doc: |-
      Additional arguments that will be passed through to bwa mem. Example value: "-k 15".
  output_file_name:
    type: string
    doc: |-
      Name of the BAM file that will be created.
arguments:
    - shellQuote: false
      valueFrom: |-
        bwa mem -t $(inputs.threads) $(inputs.args) -R "$(inputs.read_group_string)" "$(inputs.fasta_file.path)" "$(inputs.fastq_file_1.path)" "$(inputs.fastq_file_2.path)" | samtools view -b > "$(inputs.output_file_name)"
outputs:
  output_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: align_fastq_output.txt
stderr: align_fastq_error.txt
