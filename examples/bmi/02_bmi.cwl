cwlVersion: v1.1
class: CommandLineTool
doc: |-
  This tool is similar to bmi_calculator.cwl, but it retrieves the data file from a web server rather than from a local file.
requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerImageId: bmi_calculator
    dockerFile: |-
      FROM python:3.8.2-slim-buster
      RUN apt-get update && apt-get -y install wget
  NetworkAccess:
    class: NetworkAccess
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
    - entryname: calculate_bmi.py
      entry: |-
        from sys import argv
        import urllib.request

        input_file_url = argv[1]
        weight_column_name = argv[2]
        height_column_name = argv[3]
        output_file_name = argv[4]

        temp_file_path = "/tmp/biometric_data.tsv"
        urllib.request.urlretrieve(input_file_url, temp_file_path)

        with open(temp_file_path) as input_file:
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
  input_file_url:
    type: string
    doc: |-
      The URL of a tab-separated data file with biometric measurements for human individuals. The first line in the file should contain column names. One column should contain measurements representing each individual's weight (in kilograms). One column should contain measurements representing each individual's height (in centimeters). This URL must be accessible without a password.
  weight_column_name:
    type: string
    doc: |-
      Name of the column in input_file_url that contains weight measurements.
  height_column_name:
    type: string
    doc: |-
      Name of the column in input_file_url that contains height measurements.
  output_file_name:
    type: string
    doc: |-
      Name of the output file that will be created. This file will contain the same contents as input_file_url, augmented with an additional column (labeled "BMI") that indicates each person's body mass index (BMI).
arguments:
    - shellQuote: false
      valueFrom: |-
        python calculate_bmi.py "$(inputs.input_file_url)" "$(inputs.weight_column_name)" "$(inputs.height_column_name)" "$(inputs.output_file_name)"
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
stdout: 02_output.txt
stderr: 02_error.txt
