class: CommandLineTool
cwlVersion: v1.2
id: exomiser
doc: |-
  The Exomiser is a Java program that finds potential disease-causing variants
  from whole-exome or whole-genome sequencing data.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'ferlabcrsj/exomiser:2.4.1'
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
baseCommand: []
arguments:
- position: 0
  shellQuote: false
  valueFrom: >-
    tar xf $(inputs.datadir_file.path)
- position: 10
  shellQuote: false
  prefix: "&&"
  valueFrom: >-
    >&2 java -Xmx$(Math.floor(inputs.ram*1000/1.074-1))M -cp \$( cat /app/jib-classpath-file ) \$( cat /app/jib-main-class-file )
- position: 12
  shellQuote: false
  valueFrom: >-
    --output-directory=$(runtime.outdir)
- position: 15
  shellQuote: false
  valueFrom: >-
    --exomiser.data-directory=$(runtime.outdir)/data
    $(inputs.local_frequency ? ['--exomiser.', inputs.exomiser_genome, '.local-frequency-path=', inputs.local_frequency.path].join('') : "")
    $(inputs.remm_version ? ['--exomiser.remm.version="', inputs.remm_version, '" ', '--exomiser.', inputs.exomiser_genome, '.remm-path=', runtime.outdir, '/data/remm/', inputs.remm_filename].join('') : "")
    $(inputs.cadd_version ? ['--cadd.version="', inputs.cadd_version, '" ', '--exomiser.', inputs.exomiser_genome, '.cadd-snv-path=', runtime.outdir, '/data/cadd/', inputs.cadd_version, '/', inputs.cadd_snvname, ' --exomiser.', inputs.exomiser_genome, '.cadd-indel-path=', runtime.outdir, '/data/cadd/', inputs.cadd_version, '/', inputs.cadd_indelname].join('') : "")
    --exomiser.$(inputs.exomiser_genome).data-version="$(inputs.exomiser_version)"
    --exomiser.phenotype.data-version="$(inputs.exomiser_version)"
inputs:
  datadir_file: { type: "File", doc: "TAR file containing a properly formatted data directory for Exomiser." }
  vcf_file: { type: "File", secondaryFiles: [{pattern: '.tbi', required: true}], inputBinding: {prefix: "--vcf", position: 11}, doc: "Input VCF file to Exomise." }
  pheno_file: { type: "File", inputBinding: {prefix: "--sample", position: 11}, doc: "Sample file." }
  analysis_file: { type: "File", inputBinding: {prefix: "--analysis", position: 11}, doc: "Analysis file." }
  exomiser_genome: { type: "string?", default: "hg38", inputBinding: {prefix: "--assembly", position: 11}, doc: "Genome version used to generate input VCF. This value should match one of subdirectories in the data directory." }
  output_basename: { type: "string?", default: "test.exomiser", inputBinding: {prefix: "--output-filename=", separate: false, position: 12}, doc: "String to use as basename for output files." }
  output_format: { type: "string?", default: "HTML,JSON,TSV_GENE,TSV_VARIANT,VCF", inputBinding: {prefix: "--output-format=", separate: false, position: 12}, doc: "Comma-separated list of formats to output" }
  extra_args: { type: "string?", inputBinding: {position: 13}, doc: "Any extra arguments for this task." }
  exomiser_version: { type: "string", doc: "Version of exomiser" } 
  local_frequency: { type: "File?", secondaryFiles: [{pattern: '.tbi', required: true}], doc: "custom frequency source file with index." }
  remm_filename: { type: "string?", doc: "File to use for REMM annnotation. File must exist in the data_dir" }
  remm_version: { type: "string?", doc: "If using REMM, what version?" }
  cadd_snvname: { type: "string?", doc: "Filename to use for CADD SNVs. File must exist in the data_dir" }
  cadd_indelname: { type: "string?", doc: "Filename to use for CADD indels. File must exist in in the data_dir" }
  cadd_version: { type: "string?", doc: "If using CADD, what version?" }
  application_args: { type: "string?", inputBinding: {position: 16}, doc: "Any extra application arguments for this task." }
  cpu: { type: 'int?', default: 6, doc: "CPUs to allocate to this task." }
  ram: { type: 'int?', default: 36, doc: "RAM (in GB) to allocate to this task." }
outputs:
  vcf: { type: 'File?', secondaryFiles: [{pattern: '.tbi', required: true}], outputBinding: { glob: "*vcf.gz" }}
  html: { type: 'File?', outputBinding: { glob: "*html" }}
  json: { type: 'File?', outputBinding: { glob: "*json" }}
  gene_tsv: { type: 'File?', outputBinding: { glob: "*genes.tsv" }}
  variants_tsv: { type: 'File?', outputBinding: { glob: "*variants.tsv" }}
