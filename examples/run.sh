#! /bin/bash

set -o errexit

#cwl-runner --outdir Output/01_hello 01_hello.cwl hello_job.yml
#cwl-runner --outdir Output/02_hello 02_hello.cwl hello_job.yml
#cwl-runner --outdir Output/03_hello 03_hello.cwl hello_job.yml

#cwl-runner --outdir Output/BMI bmi_calculator.cwl bmi_calculator_job.yml
cwl-runner --outdir Output/BMI_url bmi_calculator_url.cwl bmi_calculator_url_job.yml
#cwl-runner --outdir Output/BMI_Ave bmi_calculator_average.cwl bmi_average_job.yml

#cwl-runner --outdir Output/DESeq2 deseq2.cwl deseq2_job.yml
