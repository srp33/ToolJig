cwlVersion: v1.2
class: CommandLineTool
label: Sum numbers from 1 file and 2 secondary files
doc: |-
  Reads a number from 1 file and 2 secondary files and then saves the sum to an output file.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: sum_files_tool
    dockerFile: |-
      FROM python:3.9-slim-buster
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
    - entryname: sum.py
      entry: |-
        import math
        import sys
        
        number_file_path1 = sys.argv[1]
        number_file_path2 = sys.argv[2]
        number_file_path3 = sys.argv[3]
        out_file_path = sys.argv[4]
        
        def readFile(file_path):
            with open(file_path) as number_file:
              return float(number_file.read().rstrip())
        
        number1 = readFile(number_file_path1)
        number2 = readFile(number_file_path2)
        number3 = readFile(number_file_path3)
        
        total = number1 + number2 + number3
        
        with open(out_file_path, 'w') as out_file:
            out_file.write(str(total))
inputs:
  number_file:
    type: File
    format: edam:format_1964
    secondaryFiles:
#      - .a
#      - .b
      - .c
      - .d
    doc: |-
      A file that contains a number in it.
  output_file:
    type: string
    doc: |-
      An output file that stores the summed value.
      #Output_File=format: edam:format_1964
arguments:
  - shellQuote: false
    valueFrom: |-
      python sum.py $(inputs.number_file.path) $(inputs.number_file.path).c $(inputs.number_file.path).d $(inputs.output_file)

      #python sum.py $(inputs.number_file.path) $(inputs.number_file.path).a $(inputs.number_file.path).b $(inputs.output_file)
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
 
s:dateCreated: "2021-03-29"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
