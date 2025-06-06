cwlVersion: v1.2
class: CommandLineTool
id: bcftools_annotate_rename_chr
doc: "Simple tool to annotate rename contigs in a VCF with a TSV old\tnew file"
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
      -o $(inputs.output_basename).$(inputs.tool_name).bcf_chr_renamed.vcf.gz
      -O z
  - position: 2
    shellQuote: false
    valueFrom: >-
      && tabix $(inputs.output_basename).$(inputs.tool_name).bcf_chr_renamed.vcf.gz
inputs:
    input_vcf: { type: 'File', secondaryFiles: [{ pattern: ".tbi", required: false }, { pattern: ".csi", required: false }],
      inputBinding: { position: 1 } }
    chr_rename_tsv: { type: 'File', doc: "tsv of old\tnew contigs",
      inputBinding: { position: 0, prefix: "--rename-chrs"} }
    threads: { type: 'int?', doc: "Number of compression/decompression threads", default: 4,
      inputBinding: { position: 0, prefix: "--threads" } }
    output_basename: string
    tool_name: string

outputs:
  bcftools_recontig_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
