cwlVersion: v1.1
class: CommandLineTool
label: Prepare a recalibration VCF file for GATK variant calling
doc: |-
  Download a recalibration VCF file from the GATK FTP server and perform some preprocessing so the VCF file will work properly, even if the reference genome was sorted differently than the VCF file.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: prep_recalibration_vcf
    dockerFile: |-
      FROM quay.io/biocontainers/picard:2.22.3--0
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
    - entryname: reconcile_vcf_with_dict.py
      entry: |-
        import sys 
        
        vcf_file_path = sys.argv[1]
        dict_file_path = sys.argv[2]
        
        with open(dict_file_path) as dict_file:
            sequences = set()
            for line in dict_file:
                if not "SN:" in line:
                    continue
        
                line_items = line.rstrip("\n").split("\t")
                sequences.add(line_items[1].replace("SN:", ""))
        
        output = ""
        with open(vcf_file_path) as vcf_file:
            for line in vcf_file:
                if line.startswith("#"):
                    if not line.startswith("##contig="):
                        output += line
                    continue
        
                if line.split("\t")[0] in sequences:
                    output += line
        
        with open(vcf_file_path, 'w') as vcf_file:
            vcf_file.write(output)
inputs:
  vcf_url:
    type: string
    doc: |-
      URL for a recalibration VCF file. Example value: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.
  fasta_file:
    type: File
    doc: |-
      Reference genome FASTA file.
    format: edam:format_1929
    secondaryFiles:
      - .dict
  output_file_name:
    type: string
    doc: |-
      Name of the output VCF file (non-gzipped). Example value: dbsnp_146.hg38.vcf.
      #Output_File=edam:format_3016
arguments:
  - shellQuote: false
    valueFrom: |-
      wget "$(inputs.vcf_url)"
      
      gunzip "`basename $(inputs.vcf_url)`"
      
      DICT_FILE="$(inputs.fasta_file.path)".dict
      
      python reconcile_vcf_with_dict.py "$(inputs.output_file_name)" "$DICT_FILE"
      
      picard -Xms128m -Xmx2g SortVcf I="$(inputs.output_file_name)" O="$(inputs.output_file_name)" SEQUENCE_DICTIONARY="$DICT_FILE" CREATE_INDEX=true
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
    format: edam:format_3016
  regular_output_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name).idx"
    doc: |-
      An index file.
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: prep_recalibration_vcf_output.txt
stderr: prep_recalibration_vcf_error.txt
 
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
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.23.owl
