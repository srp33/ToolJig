cwlVersion: v1.2
class: Workflow
id: somatic_variant_caller_workflow
label: Calls somatic variants from tumor/normal genome pair
doc: |-
  This workflow accepts FASTQ files as input and executes the steps necessary to identify somatic variants that differ between the tumor genome and the normal genome.
inputs:
  - id: download_fastq_normal_1__url
    type: string
  - id: download_fastq_normal_2__url
    type: string
  - id: download_fastq_tumor_1__url
    type: string
  - id: download_fastq_tumor_2__url
    type: string
  - id: trim_normal__args
    type: string
  - id: trim_tumor__args
    type: string
  - id: align_normal__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .fai
      - .pac
      - .sa
  - id: align_normal__threads
    type: int
  - id: align_normal__read_group_string
    type: string
  - id: align_normal__args
    type: string
  - id: align_tumor__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .fai
      - .pac
      - .sa
  - id: align_tumor__threads
    type: int
  - id: align_tumor__read_group_string
    type: string
  - id: align_tumor__args
    type: string
  - id: sort_normal__threads
    type: int
  - id: sort_tumor__threads
    type: int
  - id: mark_dups_normal__threads
    type: int
  - id: mark_dups_tumor__threads
    type: int
  - id: calculate_bqsr_normal__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
  - id: calculate_bqsr_normal__known_sites_vcf_file_1
    type: File
    format: edam:format_3016
  - id: calculate_bqsr_normal__known_sites_vcf_file_2
    type: File
    format: edam:format_3016
  - id: calculate_bqsr_normal__known_sites_vcf_file_3
    type: File
  - id: calculate_bqsr_normal__threads
    type: int
  - id: calculate_bqsr_tumor__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
  - id: calculate_bqsr_tumor__known_sites_vcf_file_1
    type: File
    format: edam:format_3016
  - id: calculate_bqsr_tumor__known_sites_vcf_file_2
    type: File
    format: edam:format_3016
  - id: calculate_bqsr_tumor__known_sites_vcf_file_3
    type: File
  - id: calculate_bqsr_tumor__threads
    type: int
  - id: apply_bqsr_normal__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
  - id: apply_bqsr_normal__threads
    type: int
  - id: apply_bqsr_tumor__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
  - id: apply_bqsr_tumor__threads
    type: int
  - id: call_small_variants__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .fai
      - .dict
  - id: call_small_variants__normal_sample_id
    type: string
  - id: call_small_variants__tumor_sample_id
    type: string
  - id: call_small_variants__threads
    type: int
  - id: call_small_variants__output_file
    type: string
  - id: call_structural_variants__fasta_file
    type: File
    format: edam:format_1929
    secondaryFiles:
      - .fai
  - id: call_structural_variants__exclude_template_url
    type: string
  - id: call_structural_variants__output_file
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
        default: temp_WwBMCvVvyj
      - id: url
        source: download_fastq_normal_1__url
    out:
      [output_file]
  download_fastq_normal_2:
    run: download_fastq_file.cwl
    in:
      - id: output_file
        default: temp_iMoHGgriPa
      - id: url
        source: download_fastq_normal_2__url
    out:
      [output_file]
  download_fastq_tumor_1:
    run: download_fastq_file.cwl
    in:
      - id: output_file
        default: temp_RVKYLHhaxR
      - id: url
        source: download_fastq_tumor_1__url
    out:
      [output_file]
  download_fastq_tumor_2:
    run: download_fastq_file.cwl
    in:
      - id: output_file
        default: temp_aFSOXnjFih
      - id: url
        source: download_fastq_tumor_2__url
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
        source: trim_normal__args
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
        source: trim_tumor__args
    out:
      [fastq_file_1, fastq_file_2]
  align_normal:
    run: align_fastq.cwl
    in:
      - id: fasta_file
        source: align_normal__fasta_file
      - id: fastq_file_1
        source: trim_normal/fastq_file_1
      - id: fastq_file_2
        source: trim_normal/fastq_file_2
      - id: threads
        source: align_normal__threads
      - id: read_group_string
        source: align_normal__read_group_string
      - id: args
        source: align_normal__args
      - id: output_file
        default: temp_kvSJywPtHt
    out:
      [output_file]
  align_tumor:
    run: align_fastq.cwl
    in:
      - id: fasta_file
        source: align_tumor__fasta_file
      - id: fastq_file_1
        source: trim_tumor/fastq_file_1
      - id: fastq_file_2
        source: trim_tumor/fastq_file_2
      - id: threads
        source: align_tumor__threads
      - id: read_group_string
        source: align_tumor__read_group_string
      - id: args
        source: align_tumor__args
      - id: output_file
        default: temp_smOvLbESau
    out:
      [output_file]
  sort_normal:
    run: sort_bam.cwl
    in:
      - id: bam_file
        source: align_normal/output_file
      - id: output_file
        default: temp_BJFdmcOkFR
      - id: threads
        source: sort_normal__threads
    out:
      [output_file]
  sort_tumor:
    run: sort_bam.cwl
    in:
      - id: bam_file
        source: align_tumor/output_file
      - id: output_file
        default: temp_SUkoDKSyvN
      - id: threads
        source: sort_tumor__threads
    out:
      [output_file]
  mark_dups_normal:
    run: mark_dups_bam.cwl
    in:
      - id: bam_file
        source: sort_normal/output_file
      - id: output_file
        default: temp_AqLXHLuSzv
      - id: threads
        source: mark_dups_normal__threads
    out:
      [output_file]
  mark_dups_tumor:
    run: mark_dups_bam.cwl
    in:
      - id: bam_file
        source: sort_tumor/output_file
      - id: output_file
        default: temp_IErhTXndAk
      - id: threads
        source: mark_dups_tumor__threads
    out:
      [output_file]
  calculate_bqsr_normal:
    run: calculate_bqsr_table.cwl
    in:
      - id: bam_file
        source: mark_dups_normal/output_file
      - id: fasta_file
        source: calculate_bqsr_normal__fasta_file
      - id: known_sites_vcf_file_1
        source: calculate_bqsr_normal__known_sites_vcf_file_1
      - id: known_sites_vcf_file_2
        source: calculate_bqsr_normal__known_sites_vcf_file_2
      - id: known_sites_vcf_file_3
        source: calculate_bqsr_normal__known_sites_vcf_file_3
      - id: output_file
        default: temp_knCJQyLxBs
      - id: threads
        source: calculate_bqsr_normal__threads
    out:
      [output_file]
  calculate_bqsr_tumor:
    run: calculate_bqsr_table.cwl
    in:
      - id: bam_file
        source: mark_dups_tumor/output_file
      - id: fasta_file
        source: calculate_bqsr_tumor__fasta_file
      - id: known_sites_vcf_file_1
        source: calculate_bqsr_tumor__known_sites_vcf_file_1
      - id: known_sites_vcf_file_2
        source: calculate_bqsr_tumor__known_sites_vcf_file_2
      - id: known_sites_vcf_file_3
        source: calculate_bqsr_tumor__known_sites_vcf_file_3
      - id: output_file
        default: temp_tNRaupRBGL
      - id: threads
        source: calculate_bqsr_tumor__threads
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
        source: apply_bqsr_normal__fasta_file
      - id: output_file
        default: temp_prdFBlNfpW
      - id: threads
        source: apply_bqsr_normal__threads
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
        source: apply_bqsr_tumor__fasta_file
      - id: output_file
        default: temp_WocnFHhuEi
      - id: threads
        source: apply_bqsr_tumor__threads
    out:
      [output_file]
  call_small_variants:
    run: call_small_variants.cwl
    in:
      - id: fasta_file
        source: call_small_variants__fasta_file
      - id: normal_bam_file
        source: apply_bqsr_normal/output_file
      - id: tumor_bam_file
        source: apply_bqsr_tumor/output_file
      - id: normal_sample_id
        source: call_small_variants__normal_sample_id
      - id: tumor_sample_id
        source: call_small_variants__tumor_sample_id
      - id: threads
        source: call_small_variants__threads
      - id: output_file
        source: call_small_variants__output_file
    out:
      [output_file]
  call_structural_variants:
    run: call_structural_variants.cwl
    in:
      - id: fasta_file
        source: call_structural_variants__fasta_file
      - id: normal_bam_file
        source: apply_bqsr_normal/output_file
      - id: tumor_bam_file
        source: apply_bqsr_tumor/output_file
      - id: exclude_template_url
        source: call_structural_variants__exclude_template_url
      - id: output_file
        source: call_structural_variants__output_file
    out:
      [output_file]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
s:dateCreated: "2021-03-03"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl