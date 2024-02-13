cwlVersion: v1.2
class: CommandLineTool
id: bcftools_annotate_vcf
doc: "Simple tool to annotate a vcf using bcftools and an annotation vcf"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest'

baseCommand: [bcftools, annotate]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      -o $(inputs.output_basename).$(inputs.tool_name).bcf_annotated.vcf.gz
      -O z
  - position: 2
    shellQuote: false
    valueFrom: >-
      && tabix $(inputs.output_basename).$(inputs.tool_name).bcf_annotated.vcf.gz
inputs:
    input_vcf: { type: 'File', secondaryFiles: ['.tbi'],
      inputBinding: { position: 1 } }
    annotation_vcf: { type: 'File', secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file",
      inputBinding: { position: 0, prefix: "--annotations"} }
    columns: { type: 'string', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF",
      inputBinding: { position: 0, prefix: "--columns" } }
    threads: { type: 'int?', doc: "Number of compression/decompression threads", default: 4,
      inputBinding: { position: 0, prefix: "--threads" } }
    output_basename: string
    tool_name: string

outputs:
  bcftools_annotated_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
