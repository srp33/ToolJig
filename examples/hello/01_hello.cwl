cwlVersion: v1.1
class: CommandLineTool
doc: |-
  This simple example tool says hello to a person by name. It uses the baseCommand property to construct the command from the inputs. (This script will only execute successfully on operating systems that provide an implementation of the "echo" command.)
baseCommand: [echo, "Hello,"]
inputs:
  given_name:
    type: string
    doc: |-
      A person's given name.
    inputBinding:
      position: 1
  surname:
    type: string
    doc: |-
      A person's surname.
    inputBinding:
      position: 2
outputs:
  standard_output:
    type: stdout
stdout: 01_output.txt
