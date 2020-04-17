#! /bin/bash

set -o errexit

#cwl-runner --outdir Output_hello1 01_hello.cwl hello_job.yml
#cwl-runner --outdir Output_hello2 02_hello.cwl hello_job.yml
#cwl-runner --outdir Output_hello3 03_hello.cwl hello_job.yml

#cwl-runner --outdir Output_BMI bmi_calculator.cwl bmi_job.yml
cwl-runner --outdir Output_BMI_Ave bmi_calculator_average.cwl bmi_average_job.yml

#cwl-runner --outdir Output_DESeq2 deseq2.cwl deseq2_job.yml
