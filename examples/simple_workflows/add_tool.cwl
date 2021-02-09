cwlVersion: v1.2
class: CommandLineTool
label: Adds two numbers
doc: |-
  This tool adds two integers and saves the sum to a file.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerImageId: add_tool
    dockerFile: |-
      FROM python:3.9-slim-buster
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
    - entryname: add.py
      entry: |-
        import sys
        
        number1 = int(sys.argv[1])
        number2 = int(sys.argv[2])
        out_file_path = sys.argv[3]
        
        with open(out_file_path, 'w') as out_file:
            total = number1 + number2
            out_file.write(str(total))
inputs:
  number1:
    type: int
    doc: |-
      An integer
  number2:
    type: int
    doc: |-
      A second integer
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created
      #Output_File=edam:format_1964
arguments:
  - shellQuote: false
    valueFrom: |-
      # This could be done with a bash command, but we're using Python for consistency with the other tool.

      python add.py $(inputs.number1) $(inputs.number2) $(inputs.output_file_name)
outputs:
  output_from_input_1:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Output file matching the name specified in the "output_file_name" input.
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
