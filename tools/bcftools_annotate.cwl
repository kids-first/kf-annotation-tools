cwlVersion: v1.2
class: CommandLineTool
id: bcftools_annotate_vcf
doc: >-
  Simple tool to modify annotations in a VCF by doing one or more of the following:
    - Add annotations from another VCF based on position matching
      - Can specify which columns to add from the annotation VCF
    - Strip existing INFO annotations
    - Rename contigs based on a provided TSV of old\tnew names
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/bcftools:1.20'

baseCommand: [bcftools, annotate]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      -o $(inputs.output_basename).$(inputs.tool_name).bcf_annotated.vcf.gz
      -O z
      --write-index=tbi

inputs:
    input_vcf: { type: 'File', secondaryFiles: ['.tbi'],
      inputBinding: { position: 1 } }
    annotation_vcf: { type: 'File?', secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file",
      inputBinding: { position: 0, prefix: "--annotations"} }
    columns: { type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF",
      inputBinding: { position: 0, prefix: "--columns" } }
    threads: { type: 'int?', doc: "Number of compression/decompression threads", default: 4,
      inputBinding: { position: 0, prefix: "--threads" } }
    strip_info: {type: 'string?', doc: "If given, remove previous annotation information based on INFO field, i.e. to strip VEP info, use INFO/ANN",
      inputBinding: { position: 0, prefix: "-x" } }
    chr_rename_tsv: { type: 'File?', doc: "tsv of old\tnew contigs",
      inputBinding: { position: 0, prefix: "--rename-chrs"} }
    output_basename: string
    tool_name: string

outputs:
  bcftools_annotated_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
