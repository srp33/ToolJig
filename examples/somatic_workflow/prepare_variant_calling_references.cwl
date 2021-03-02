cwlVersion: v1.2
class: Workflow
id: prepare_variant_calling_references
label: Prepare reference files for somatic variant calling
doc: |-
  This workflow downloads and prepares reference-genome and recalibration files to be used in somatic variant calling.
inputs:
  - id: prep_ref_genome__ref_genome_version
    type: string
  - id: prep_recalibration_vcf_1000G__output_file
    type: string
  - id: prep_recalibration_vcf_1000G__vcf_url
    type: string
  - id: prep_recalibration_vcf_dbsnp__output_file
    type: string
  - id: prep_recalibration_vcf_dbsnp__vcf_url
    type: string
  - id: prep_recalibration_vcf_indels__output_file
    type: string
  - id: prep_recalibration_vcf_indels__vcf_url
    type: string
outputs:
  - id: prep_recalibration_vcf_1000G__output_file
    type: File
    outputSource: prep_recalibration_vcf_1000G/output_file
  - id: prep_recalibration_vcf_dbsnp__output_file
    type: File
    outputSource: prep_recalibration_vcf_dbsnp/output_file
  - id: prep_recalibration_vcf_indels__output_file
    type: File
    outputSource: prep_recalibration_vcf_indels/output_file
steps:
  prep_ref_genome:
    run: prep_ref_genome.cwl
    in:
      - id: ref_genome_version
        source: prep_ref_genome__ref_genome_version
    out:
      [reference_fasta_file]
  prep_recalibration_vcf_1000G:
    run: prep_recalibration_vcf.cwl
    in:
      - id: fasta_file
        source: prep_ref_genome/reference_fasta_file
      - id: output_file
        source: prep_recalibration_vcf_1000G__output_file
      - id: vcf_url
        source: prep_recalibration_vcf_1000G__vcf_url
    out:
      [output_file]
  prep_recalibration_vcf_dbsnp:
    run: prep_recalibration_vcf.cwl
    in:
      - id: fasta_file
        source: prep_ref_genome/reference_fasta_file
      - id: output_file
        source: prep_recalibration_vcf_dbsnp__output_file
      - id: vcf_url
        source: prep_recalibration_vcf_dbsnp__vcf_url
    out:
      [output_file]
  prep_recalibration_vcf_indels:
    run: prep_recalibration_vcf.cwl
    in:
      - id: fasta_file
        source: prep_ref_genome/reference_fasta_file
      - id: output_file
        source: prep_recalibration_vcf_indels__output_file
      - id: vcf_url
        source: prep_recalibration_vcf_indels__vcf_url
    out:
      [output_file]
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
s:dateCreated: "2021-03-02"
s:license: https://spdx.org/licenses/Apache-2.0
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl