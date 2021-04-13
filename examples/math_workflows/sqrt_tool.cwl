cwlVersion: v1.2
class: CommandLineTool
label: Calculates the square root of a number
doc: |-
  This tool reads a number from a file, calculates the square root of that number, and saves the output to a file.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: sqrt_tool
    dockerFile: |-
      FROM python:3.9-slim-buster
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
    - entryname: sqrt.py
      entry: |-
        import math
        import sys
        
        number_file_path = sys.argv[1]
        out_file_path = sys.argv[2]
        
        with open(number_file_path) as file1:
            number = float(file1.read())
        
        with open(out_file_path, 'w') as out_file:
            out_file.write(str(math.sqrt(number)))
inputs:
  number_file:
    type: File
    doc: |-
      A file with a single number in it
    format: edam:format_1964
  output_file:
    type: string
    doc: |-
      Name of the output file
      #Output_File=edam:format_1964
arguments:
  - shellQuote: false
    valueFrom: |-
      python sqrt.py $(inputs.number_file.path) $(inputs.output_file)
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file)"
    doc: |-
      Output file matching the name specified in the "output_file" input.
    format: edam:format_1964
  standard_output:
    type: stdout
    format: edam:format_1964
  standard_error:
    type: stderr
    format: edam:format_1964
stdout: stdout.txt
stderr: stderr.txt
 
s:author:
  - class: s:Person
    s:name: Stephen Piccolo
    s:identifier: https://orcid.org/0000-0003-2001-5640
 
s:dateCreated: "2020-07-14"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
