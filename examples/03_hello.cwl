cwlVersion: v1.1
class: CommandLineTool
doc: |-
  This simple example tool says hello to a person by name. This example is identical to hello2.cwl, except that the tool is executed within a container environment. We use a container image pulled from Docker Hub that consists of the "buster" release (version 10.3) of the Debian Linux operating system. It is a "slim" version of the operating system, meaning that it provides only essential components.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: debian:buster-slim
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
        echo Hello, $(inputs.given_name) $(inputs.surname)! You are $(inputs.age) years old.
outputs:
  standard_output:
    type: stdout
stdout: output.txt
