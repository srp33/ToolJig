cwlVersion: v1.1
class: CommandLineTool
doc: |
  This tool uses an embedded Python script to find the body mass index (BMI) of humans whose weights and heights are stored in a tab-separated data file. It then calculates and prints the average BMI for these individuals. This example demonstrates a way to support reproducibility of computational analyses: we embed within the tool the data file that is used for the analysis.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: python:3.8.2-slim-buster
  InitialWorkDirRequirement:
    listing:
    - entryname: calculate_bmi_average.py
      entry: |-
        from sys import argv

        output_file_name = argv[1]

        bmi_values = []

        with open("biometric_data.tsv") as input_file:
            header_items = input_file.readline().rstrip("\n").split("\t")
            weight_i = header_items.index("Weight")
            height_i = header_items.index("Height")

            for line in input_file:
                line_items = line.rstrip("\n").split("\t")
                weight_kg = float(line_items[weight_i])
                height_cm = float(line_items[height_i])
                bmi_values.append((weight_kg / height_cm**2) * 10000)

        average_bmi = sum(bmi_values) / len(bmi_values)

        with open(output_file_name, "w") as output_file:
            output_file.write(f"{average_bmi:.2f}")

        print("Calculations were completed successfully!")
    - entryname: biometric_data.tsv
      entry: |-
        Surname	GivenName	Weight	Height
        Carson	Hafsa	82.4	167.4
        Alvarado	Yasmine	69.2	172.4
        Green	Vanessa	66.8	194.0
inputs:
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created. This file will contain the average BMI value, rounded to two decimal places.
arguments:
    - shellQuote: false
      valueFrom: |-
        python calculate_bmi_average.py "$(inputs.output_file_name)"
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file_name)"
    doc: |-
      Here we indicate that an output file matching the name specified in the inputs should be generated.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: output.txt
stderr: error.txt
