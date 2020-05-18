cwlVersion: v1.1
class: CommandLineTool
doc: |-
  This simple example tool says hello to a person by name. It uses the arguments property to construct the command from the inputs, thus providing more flexibility for constructing the command (we can easily put an exclamation point at the end of the greeting). Note also that we do not need to specify "position" input bindings. (This script will only execute successfully on operating systems that provide an implementation of the "echo" command.)
requirements:
  ShellCommandRequirement: {}
inputs:
  given_name:
    type: string
    doc: |-
      A person's given name.
  surname:
    type: string
    doc: |-
      A person's surname.
  age:
    type: int
    doc: |-
      A person's age in years.
arguments:
    - shellQuote: false
      valueFrom: |-
        echo "Hello, $(inputs.given_name) $(inputs.surname)! You are $(inputs.age) years old."
outputs:
  standard_output:
    type: stdout
stdout: 02_output.txt
