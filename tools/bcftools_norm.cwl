cwlVersion: v1.2
class: CommandLineTool
id: bcftools_norm
doc: "Normalize VCF for DB loading and VCF-to-VCF variant comparison"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/bcftools:1.20'
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  - class: InlineJavascriptRequirement

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      bcftools norm
      --write-index=tbi
      --old-rec-tag OLD_RECORD
      -c w
      -m -any
      -o $(inputs.output_basename).$(tool_name).vcf.gz

inputs:
  input_vcf: { type: File, inputBinding: { position: 10 } }
  threads: { type: 'int?', default: 4, inputBinding: { position: 1, prefix: "--threads" } }
  output_type: { type: [ 'null', {type: enum, name: output_type, symbols: [ "u", "b", "v", "z"]}], default: "z", inputBinding: { position: 1, prefix: "-O" } }
  fasta: { type: File, inputBinding: { position: 1, prefix: "-f" } }
  output_basename: string
  tool_name: string
outputs:
  normalized_vcf:
    type: File
    outputBinding:
      glob: "*.{v,b}cf{,.gz}"
    secondaryFiles: ['.tbi?']
