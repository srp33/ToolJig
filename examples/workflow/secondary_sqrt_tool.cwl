cwlVersion: v1.2
class: CommandLineTool
label: Calculates the square root of a number plus numbers stored in secondary files
doc: |-
  This tool reads an integer from a file, calculates the square root of that number, and saves it to an output file. It does the same for two secondary files.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: secondary_sqrt_tool
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
        
        with open(number_file_path) as number_file:
            with open(out_file_path, 'w') as out_file:
              out_file.write(str(math.sqrt(float(number_file.read()))))

        with open(number_file_path + ".a") as number_file:
            with open(out_file_path + ".c", 'w') as out_file:
              out_file.write(str(math.sqrt(float(number_file.read()))))

        with open(number_file_path + ".b") as number_file:
            with open(out_file_path + ".d", 'w') as out_file:
              out_file.write(str(math.sqrt(float(number_file.read()))))

inputs:
  number_file:
    type: File
    doc: |-
      A file with a single number in it
    format: edam:format_1964
    secondaryFiles:
      - .a
      - .b
  output_file_name:
    type: string
    doc: |-
      Name of the output file
      #Output_File=edam:format_1964
arguments:
  - shellQuote: false
    valueFrom: |-
      python sqrt.py $(inputs.number_file.path) $(inputs.output_file_name)
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
    format: edam:format_1964
    secondaryFiles:
      - .c
      - .d
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
