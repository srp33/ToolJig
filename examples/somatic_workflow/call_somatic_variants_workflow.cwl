cwlVersion: v1.2
class: Workflow
id: call_somatic_variants_workflow
label: Calls somatic variants from tumor/normal genome pair
doc: |-
  This workflow accepts (gzipped) FASTQ files as input and executes the steps necessary to identify somatic variants that differ between the tumor genome and the normal genome.
inputs:
  - id: fastq_normal_1_url
    type: string
  - id: fastq_normal_2_url
    type: string
  - id: fastq_tumor_1_url
    type: string
  - id: fastq_tumor_2_url
    type: string
  - id: trim_args
    type: string
  - id: fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .fai
      - .pac
      - .sa
  - id: threads
    type: int
  - id: read_group_string_normal
    type: string
  - id: align_args
    type: string
  - id: read_group_string_tumor
    type: string
  - id: known_sites_vcf_file_1
    type: File
    format: edam:format_3016
  - id: known_sites_vcf_file_2
    type: File
    format: edam:format_3016
  - id: known_sites_vcf_file_3
    type: File
    format: edam:format_3016
  - id: normal_sample_id
    type: string
  - id: tumor_sample_id
    type: string
  - id: small_variants_output_file
    type: string
  - id: delly_exclude_template_url
    type: string
  - id: sv_output_file
    type: string
outputs:
  - id: call_small_variants__output_file
    type: File
    outputSource: call_small_variants/output_file
  - id: call_structural_variants__output_file
    type: File
    outputSource: call_structural_variants/output_file
steps:
  download_fastq_normal_1:
    run: download_fastq_file.cwl
    in:
      - id: output_file
        default: fastq_normal_1.fastq.gz
      - id: url
        source: fastq_normal_1_url
    out:
      [output_file]
  download_fastq_normal_2:
    run: download_fastq_file.cwl
    in:
      - id: output_file
        default: fastq_normal_2.fastq.gz
      - id: url
        source: fastq_normal_2_url
    out:
      [output_file]
  download_fastq_tumor_1:
    run: download_fastq_file.cwl
    in:
      - id: output_file
        default: fastq_tumor_1.fastq.gz
      - id: url
        source: fastq_tumor_1_url
    out:
      [output_file]
  download_fastq_tumor_2:
    run: download_fastq_file.cwl
    in:
      - id: output_file
        default: fastq_tumor_2.fastq.gz
      - id: url
        source: fastq_tumor_2_url
    out:
      [output_file]
  trim_normal:
    run: trim_fastq.cwl
    in:
      - id: fastq_file_1
        source: download_fastq_normal_1/output_file
      - id: fastq_file_2
        source: download_fastq_normal_2/output_file
      - id: args
        source: trim_args
    out:
      [fastq_file_1, fastq_file_2]
  trim_tumor:
    run: trim_fastq.cwl
    in:
      - id: fastq_file_1
        source: download_fastq_tumor_1/output_file
      - id: fastq_file_2
        source: download_fastq_tumor_2/output_file
      - id: args
        source: trim_args
    out:
      [fastq_file_1, fastq_file_2]
  align_normal:
    run: align_fastq.cwl
    in:
      - id: fasta_file
        source: fasta_file
      - id: fastq_file_1
        source: trim_normal/fastq_file_1
      - id: fastq_file_2
        source: trim_normal/fastq_file_2
      - id: threads
        source: threads
      - id: read_group_string
        source: read_group_string_normal
      - id: args
        source: align_args
      - id: output_file
        default: aligned_normal.bam
    out:
      [output_file]
  align_tumor:
    run: align_fastq.cwl
    in:
      - id: fasta_file
        source: fasta_file
      - id: fastq_file_1
        source: trim_tumor/fastq_file_1
      - id: fastq_file_2
        source: trim_tumor/fastq_file_2
      - id: threads
        source: threads
      - id: read_group_string
        source: read_group_string_tumor
      - id: args
        source: align_args
      - id: output_file
        default: aligned_tumor.bam
    out:
      [output_file]
  sort_normal:
    run: sort_bam.cwl
    in:
      - id: bam_file
        source: align_normal/output_file
      - id: output_file
        default: sorted_normal.bam
      - id: threads
        source: threads
    out:
      [output_file]
  sort_tumor:
    run: sort_bam.cwl
    in:
      - id: bam_file
        source: align_tumor/output_file
      - id: output_file
        default: sorted_tumor.bam
      - id: threads
        source: threads
    out:
      [output_file]
  mark_dups_normal:
    run: mark_dups_bam.cwl
    in:
      - id: bam_file
        source: sort_normal/output_file
      - id: output_file
        default: dups_marked_normal.bam
      - id: threads
        source: threads
    out:
      [output_file]
  mark_dups_tumor:
    run: mark_dups_bam.cwl
    in:
      - id: bam_file
        source: sort_tumor/output_file
      - id: output_file
        default: dups_marked_tumor.bam
      - id: threads
        source: threads
    out:
      [output_file]
  calculate_bqsr_normal:
    run: calculate_bqsr_table.cwl
    in:
      - id: bam_file
        source: mark_dups_normal/output_file
      - id: fasta_file
        source: fasta_file
      - id: known_sites_vcf_file_1
        source: known_sites_vcf_file_1
      - id: known_sites_vcf_file_2
        source: known_sites_vcf_file_2
      - id: known_sites_vcf_file_3
        source: known_sites_vcf_file_3
      - id: output_file
        default: random__iVBcjXmFTy
      - id: threads
        source: threads
    out:
      [output_file]
  calculate_bqsr_tumor:
    run: calculate_bqsr_table.cwl
    in:
      - id: bam_file
        source: mark_dups_tumor/output_file
      - id: fasta_file
        source: fasta_file
      - id: known_sites_vcf_file_1
        source: known_sites_vcf_file_1
      - id: known_sites_vcf_file_2
        source: known_sites_vcf_file_2
      - id: known_sites_vcf_file_3
        source: known_sites_vcf_file_3
      - id: output_file
        default: random__TNHumTCAVr
      - id: threads
        source: threads
    out:
      [output_file]
  apply_bqsr_normal:
    run: apply_bqsr_bam.cwl
    in:
      - id: bam_file
        source: mark_dups_normal/output_file
      - id: bqsr_table_file
        source: calculate_bqsr_normal/output_file
      - id: fasta_file
        source: fasta_file
      - id: output_file
        default: final_normal.bam
      - id: threads
        source: threads
    out:
      [output_file]
  apply_bqsr_tumor:
    run: apply_bqsr_bam.cwl
    in:
      - id: bam_file
        source: mark_dups_tumor/output_file
      - id: bqsr_table_file
        source: calculate_bqsr_tumor/output_file
      - id: fasta_file
        source: fasta_file
      - id: output_file
        default: final_tumor.bam
      - id: threads
        source: threads
    out:
      [output_file]
  call_small_variants:
    run: call_small_variants.cwl
    in:
      - id: fasta_file
        source: fasta_file
      - id: normal_bam_file
        source: apply_bqsr_normal/output_file
      - id: tumor_bam_file
        source: apply_bqsr_tumor/output_file
      - id: normal_sample_id
        source: normal_sample_id
      - id: tumor_sample_id
        source: tumor_sample_id
      - id: threads
        source: threads
      - id: output_file
        source: small_variants_output_file
    out:
      [output_file]
  call_structural_variants:
    run: call_structural_variants.cwl
    in:
      - id: fasta_file
        source: fasta_file
      - id: normal_bam_file
        source: apply_bqsr_normal/output_file
      - id: tumor_bam_file
        source: apply_bqsr_tumor/output_file
      - id: exclude_template_url
        source: delly_exclude_template_url
      - id: output_file
        source: sv_output_file
    out:
      [output_file]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
s:dateCreated: "2021-04-08"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl