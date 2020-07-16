cwlVersion: v1.1
class: CommandLineTool
label: add tool
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: add_tool
    dockerFile: |-
      FROM python3.8.2
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
      A second number
  output_file_name:
    type: string
    doc: |-
      Name of the output file
      #Output_File=edam:format_1964
arguments:
  - shellQuote: false
    valueFrom: |-
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
 
s:dateCreated: "2020-07-14"
s:license: https://spdx.org/licenses/Apache-2.0
 
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.23.owl
