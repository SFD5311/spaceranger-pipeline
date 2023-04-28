cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/visium-quantification:spaceranger
  ResourceRequirement:
    ramMin: 28672
baseCommand: /opt/space_ranger_entrypoint.py
label: Run spaceranger tool on FASTQ input

# arguments are hardcoded in salmon_wrapper.py

inputs:
  fastq_dir:
    type: Directory
    inputBinding:
      position: 1
  img_path:
    type: File
    inputBinding:
      position: 2
  gpr_path:
    type: File
    inputBinding:
      position: 3
  assay:
    type: string
    inputBinding:
      position: 4

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: outs
