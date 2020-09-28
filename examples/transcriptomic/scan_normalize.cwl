cwlVersion: v1.1
class: CommandLineTool
label: SCAN normalization example
doc: |-
  This tool demonstrates how to normalize gene-expression microarray data stored in Gene Expression Omnibus and create a graph illustrating technical artifacts. It uses the Bioconductor SCAN.UPC package for normalization. We use a container image that provides core Bioconductor components (release version 3.10) and use R code to install this package.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: scan_normalize
    dockerFile: |-
      FROM bioconductor/bioconductor_docker:RELEASE_3_11
      
      RUN R -e 'BiocManager::install("SCAN.UPC")'
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
    - entryname: normalize.R
      entry: |-
        library(SCAN.UPC)
        
        # Parse the command-line arguments
        geo_series_id = commandArgs()[8]
        out_file_name = commandArgs()[9]
        
        # Download the data, perform normalization, and save the data to a file
        data = SCAN(geo_series_id, outFilePath=out_file_name)
inputs:
  geo_series_id:
    type: string
    doc: |-
      The unique identifier for a GEO series with Affymetrix microarray data.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created. This will be a tab-separated file.
      #Output_File=edam:format_3751
arguments:
  - shellQuote: false
    valueFrom: |-
      export R_LIBS="/tmp"
      
      Rscript normalize.R "$(inputs.geo_series_id)" "$(inputs.output_file_name)"
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
    format: edam:format_3751
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: scan_output.txt
stderr: scan_error.txt
 
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
