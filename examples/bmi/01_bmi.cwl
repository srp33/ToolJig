cwlVersion: v1.2
class: CommandLineTool
doc: |-
  This tool demonstrates a way to process an input data file using an embedded Python script. The script requires as input a data file with tab-separated values. Two of the columns in this file must represent weights and heights of human individuals. The script calculates the individuals' body mass index (BMI) values, stores them in a new column in the file, and writes the file to the output directory. We use a container image from Docker Hub that is based on the "buster" release (version 10.3) of the Debian Linux operating system and includes version 3.8.2 of the Python interpreter.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: python:3.9-slim-buster
  InitialWorkDirRequirement:
    listing:
    - entryname: calculate_bmi.py
      entry: |-
        from sys import argv

        input_file_path = argv[1]
        weight_column_name = argv[2]
        height_column_name = argv[3]
        output_file_name = argv[4]

        with open(input_file_path) as input_file:
            with open(output_file_name, "w") as output_file:
                header_items = input_file.readline().rstrip("\n").split("\t")
                weight_i = header_items.index(weight_column_name)
                height_i = header_items.index(height_column_name)

                output_file.write("\t".join(header_items + ["BMI"]) + "\n")

                for line in input_file:
                    line_items = line.rstrip("\n").split("\t")
                    weight_kg = float(line_items[weight_i])
                    height_cm = float(line_items[height_i])
                    bmi = (weight_kg / height_cm**2) * 10000

                    line_items.append(f"{bmi:.1f}")
                    output_file.write("\t".join(line_items) + "\n")

        print("Calculations were completed successfully!")
inputs:
  input_file:
    type: File
    doc: |-
      A tab-separated data file with biometric measurements for human individuals. The first line in the file should contain column names. One column should contain measurements representing each individual's weight (in kilograms). One column should contain measurements representing each individual's height (in centimeters).
  weight_column_name:
    type: string
    doc: |-
      Name of the column in input_file that contains weight measurements.
  height_column_name:
    type: string
    doc: |-
      Name of the column in input_file that contains height measurements.
  output_file:
    type: string
    doc: |-
      Name of the output file that will be created. This file will contain the same contents as input_file, augmented with an additional column (labeled "BMI") that indicates each person's body mass index (BMI).
arguments:
    - shellQuote: false
      valueFrom: |-
        python calculate_bmi.py "$(inputs.input_file.path)" "$(inputs.weight_column_name)" "$(inputs.height_column_name)" "$(inputs.output_file)"
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_file)"
    doc: |-
      Output file matching the name specified in the "output_file" input.
  standard_output:
    type: stdout
  standard_error:
    type: stderr
stdout: 01_output.txt
stderr: 01_error.txt
