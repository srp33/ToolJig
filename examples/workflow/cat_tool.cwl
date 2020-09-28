cwlVersion: v1.1
class: CommandLineTool
label: cat tool
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: cat_tool
    dockerFile: |-
      FROM python3.8.2
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
    - entryname: cat.py
      entry: |-
        import sys
        
        file_path1 = sys.argv[1]
        file_path2 = sys.argv[2]
        out_file_path = sys.argv[3]
        
        with open(file_path1) as file1:
            contents1 = file1.read()
        
        with open(file_path2) as file2:
            contents2 = file2.read()
        
        with open(out_file_path, 'w') as out_file:
            out_file.write(contents1 + "\n")
            out_file.write(contents2 + "\n")
inputs:
  file1:
    type: File
    doc: |-
      An file with text
    format: edam:format_1964
  file2:
    type: File
    doc: |-
      A second file
    format: edam:format_1964
  output_file_name:
    type: string
    doc: |-
      Name of the output file
      #Output_File=edam:format_1964
arguments:
  - shellQuote: false
    valueFrom: |-
      python cat.py $(inputs.file1.path) $(inputs.file2.path) $(inputs.output_file_name)
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
 - https://schema.org/version/latest/schemaorg-current-http.rdf
 - http://edamontology.org/EDAM_1.23.owl
