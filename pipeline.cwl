#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: scRNA-seq pipeline using Salmon and Alevin
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
inputs:
  fastq_dir:
    label: "Directory containing FASTQ files"
    type: Directory
  img_path:
    label: "Directory containing TIFF files"
    type: Directory
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
outputs:
  output_directory:
    outputSource: space-ranger/output_directory
    type: Directory
    label: "Directory containing spaceranger outputs"

steps:
  space-ranger:
    in:
      fastq_dir:
        source: fastq_dir
      img_path:
        source: img_path

      threads:
        source: threads
    out:
      - h5ad_file
      - bam_file
    run: steps/quantification.cwl
