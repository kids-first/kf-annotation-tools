cwlVersion: v1.2
class: Workflow
id: kfdrc-somatic-snv-annot-wf
label: Kids First DRC Somatic SNV Annotation Workflow
doc: |
  # Kids First DRC Somatic Variant Annotation Workflow
  This is a subworkflow that is part of the Kids First DRC Somatic Variant Workflow that can be run as standalone.
  Annotation of variant calls helps give context to the possible biological consequences of each variant.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  It does the following things as described below:

  1. Normalize VCF
  1. Strip specified `INFO` and `FORMAT` fields (Only if adding a new annotation that clashes with existing)
  1. Annotate with VEP - can be skipped if VEP run previously and downstream tools are to be repeated
  1. Annotated with an additional vcf - optional, recommend using a gnomAD VCF with at least AF
  1. Soft filter on remarkable variant characteristics
     - KF recommends normal read depth <= 7 and gnomAD AF > 0.001 and gnomAD FILTER == PASS
     - This output will be considered `protected`
  1. Annotate with hotspots - KF recommends cancer genome hotspots v2, formatting required and explained below
  1. Create MAF output using a modified version of MSKCC's vcf2maf
  1. Hard filter on vcf based on user-specified criteria - this output would be considered `public`
  1. Create MAF output based on `public` vcf

  ![annot workflow flowchart](../docs/somatic_annotation_wf.png)

  ## Workflow Description and KF Recommended Inputs
  The additional gnomAD annotation, hotspot annotation, and soft + hard filtering are part of process called "Germline Masking."
  The purpose of this is to create outputs that are safe for public consumption by marking and creating a version of outputs deemed a "germline risk" based on specified criteria.
  For KF, based on guidance from the Genomic Data Commons (GDC), this means filtering variants with a normal read depth of <= 7 reads and an AF score greater than 0.001 from a PASS gnomAD record.
  The gnomAD AF filter is pretty intuitive - gnomAD is a database resource of variants and their estimated prevalence in the human population.
  Therefore, a variant that is higher than the recommended threshold can be seen as a higher risk of being a common and identifiable variant, and a lower risk for being disease-causing.
  Additionally, this AF score needs to come from a PASS/non-filtered record. Therefore, records that do not have their gnomad_FILTER set to PASS are not considered.
  The normal depth argument may be less intuitive to some, and an example might help explain its importance:
  You've got a somatic variant candidate, and the normal-sample reads all support the reference allele.
  However,  if there are only ~5 of those reads, then there's a decent chance that it's a germline het variant and you just didn't get lucky enough to see any alt reads in the normal sample.
  A prior understanding that heterozygous germline variant are much more common than somatic variants informs this.

  ### Recommended reference inputs - all file references can be obtained [here](https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/)
  Secondary files needed for each reference file will be a sub-bullet point
   - `indexed_reference_fasta`: `Homo_sapiens_assembly38.fasta`
     - `Homo_sapiens_assembly38.fasta.fai`
     - `Homo_sapiens_assembly38.dict`
   - `echtvar_anno_zips`: `gnomad.v3.1.1.custom.echtvar.zip`
   - `bcftools_strip_columns`: csv string of columns to strip if needed to avoid conflict, i.e INFO/AF
   - `bcftools_public_filter`: 
     - DGD nexus export: `FILTER="OK;clinicalReported"|FILTER="clinicalReported"`
     - All others: `FILTER="PASS"|INFO/HotSpotAllele=1`
   - `gatk_filter_name`:
    - DGD nexus export: null
    - All others: ["NORM_DP_LOW", "GNOMAD_AF_HIGH"]
   - `gatk_filter_expression`:
     - DGD nexus export: null
     - All others: [`vc.getGenotype('insert_normal_sample_name').getDP() <= 7`, `gnomad_3_1_1_AF != '.' && gnomad_3_1_1_AF > 0.001 && gnomad_3_1_1_FILTER=='PASS'`] # NOTE!! Replace `insert_normal_sample_name` with the value you'd use for `input_normal_name`! # NOTE!! If your annotation includes dot values, those values must first be excluded! If they are not, GATK will error trying to convert those values!
   - `vep_cache`: `homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz`
   - `genomic_hotspots`: `tert.bed` # This file has two common TERT promoter gene hot spots
   - `protein_snv_hotspots`: `kfdrc_protein_snv_cancer_hotspots_20240718.txt` #  Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots. File header contains generation history
   - `protein_indel_hotspots`: `protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv` # A tsv formatted INDEL subset of https://www.cancerhotspots.org/files/hotspots_v2.xls
   - `custom_enst`: `kf_isoform_override.tsv` # As of VEP 104, several genes have had their canonical transcripts redefined. While the VCF will have all possible isoforms, this affects maf file output and may results in representative protein changes that defy historical expectations

  ### Source-specific inputs
  For each input, the sub-bullet refers to when to use the suggested input
   - `add_common_fields`
     - Strelka2 calls: `true`, *exception if already run previously and other downstream tools are being run*
     - All others: `false`
   - `bcftools_recontig_tsv`: _DGD nexus export ONLY_: For inputs with chr stripped, provide TSV with `old\tnew` contigs
   - `bcftools_prefilter_csv`: _DGD nexus export ONLY_: `FILTER="OK;clinicalReported"|FILTER="OK"|FILTER="clinicalReported"`
   - `retain_info` # This is fairly subjective, some useful columns unique from each caller to carry over from VCF to MAF
     - Strelka2:
        ```
        gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MQ,MQ0,QSI,HotSpotAllele
        ```
     - Mutect2:
        ```
        gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MBQ,TLOD,HotSpotAllele
        ```
     - Lancet:
        ```
        gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MS,FETS,HotSpotAllele
        ```
     - Vardict:
        ```
        gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MSI,MSILEN,SOR,SSF,HotSpotAllele
        ```
     - Consensus:
        ```
        gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MQ,MQ0,CAL,HotSpotAllele
        ```
     - DGD Nexus:
       ```
       gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,Classification,GenomicSource,ClinicallyReported,ManuallyEntered,Correlation,HotSpotAllele
       ```
   - `retain_ann` # Similar to above, if run for KF harmonization, recommend the following:
     - ALL: `HGVSg`
  - `bcftools_strip_columns` # if reannotating an old file:
     - `FILTER/GNOMAD_AF_HIGH,FILTER/NORM_DP_LOW,INFO/CSQ,INFO/HotSpotAllele` # recommended if re-annotating from an older VEP cache
     - `FILTER/GNOMAD_AF_HIGH,FILTER/NORM_DP_LOW,INFO/HotSpotAllele` # recommended if repeating hot spot and want to keep VEP
   - `bcftools_prefilter_csv` # if annotating a file with calls you want screen for, use this. i.e `FILTER="PASS"`
   - `disable_norm` # set to `True` if existing input already normalized or of you have justification for skipping this step
   - `disable_vep_annotation` # set to `True` if existing VEP annotation of file is ok
   - `disable_hotspot_annotation` # set to `True` if existing HotSpot annotation is ok
   - `tool_name`:
     - `Strelka2`: `strelka2_somatic`
     - `Mutect2`: `mutect2_somatic`
     - `Lancet`: `lancet_somatic`
     - `VarDict Java`: `vardict_somatic`
     - `consensus`: `consensus_somatic`
     - `DGD nexus export`: `dgd_nexus`
   - `vep_cores`: 
     - DGD nexus export: `8`
     - Otherwise: `16`
   - `vep_ram`: 
     - DGD nexus export: `8`
     - Otherwise: `32`
   - `vep_buffer`: `5000`

  ## Workflow outputs
   - `annotated_protected`: `PASS` VCF with annotation pipeline soft `FILTER`-added values, VCF index, and MAF format of VCF
   - `annotated_public_vcf`: Same as `annotated_protected`, hard-filtered to include `PASS` only
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict], "sbg:suggestedValue": {class: File, path: 60639014357c3a53540ca7a3,
      name: Homo_sapiens_assembly38.fasta, secondaryFiles: [{class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai},
        {class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict}]}}
  input_vcf: {type: 'File', secondaryFiles: [{pattern: ".tbi", required: false}, {pattern: ".csi", required: false}], doc: "Input
      vcf to annotate and soft filter"}
  input_tumor_name: string
  input_normal_name: string
  add_common_fields: {type: 'boolean', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  bcftools_prefilter_csv: {type: 'string?', doc: "csv of bcftools filter params if you want to prefilter before annotation"}
  bcftools_recontig_tsv: {type: 'File?', doc: "TSV file of old\tnew contigs, if needed"}
  disable_norm: {type: 'boolean?', doc: "Skip normalization step. Not recommended unless input is already normalized", default: false}
  disable_hotspot_annotation: {type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task."}
  disable_vep_annotation: {type: 'boolean?', doc: "Disable VEP Annotation and skip this task.", default: false}
  echtvar_anno_zips: {type: 'File[]?', doc: "Annotation ZIP files for echtvar anno", "sbg:suggestedValue": [{class: File, path: 65c64d847dab7758206248c6,
        name: gnomad.v3.1.1.custom.echtvar.zip}]}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]?', doc: "Array of names for each filter tag to add"}
  gatk_filter_expression: {type: 'string[]?', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration
      for clues"}
  vep_ram: {type: 'int?', default: 32, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', default: 16, doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', default: 1000, doc: "Increase or decrease to balance speed and memory usage"}
  vep_cache: {type: 'File?', doc: "tar gzipped cache from ensembl/local converted cache", "sbg:suggestedValue": {class: File, path: 6332f8e47535110eb79c794f,
      name: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz}}
  dbnsfp: {type: 'File?', secondaryFiles: [.tbi, ^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing
      dbNSFP annotations"}
  dbnsfp_fields: {type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all"}
  merged: {type: 'boolean?', doc: "Set to true if merged cache used", default: true}
  run_cache_existing: {type: 'boolean?', doc: "Run the check_existing flag for cache", default: true}
  run_cache_af: {type: 'boolean?', doc: "Run the allele frequency flags for cache", default: true}
  run_stats: {type: 'boolean?', doc: "Create stats file? Disable for speed", default: false}
  cadd_indels: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations"}
  cadd_snvs: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations"}
  genomic_hotspots: {type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to
      hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a423, name: tert.bed}]}
  protein_snv_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid
      positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 66980e845a58091951d53984, name: kfdrc_protein_snv_cancer_hotspots_20240718.txt}]}
  protein_indel_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino
      acid position ranges corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 663d2bcc27374715fccd8c6f, name: protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv}]}
  output_basename: string
  tool_name: string
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  retain_ann: {type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF"}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}
  custom_enst: {type: 'File?', doc: "Use a file with ens tx IDs for each gene to override VEP PICK", "sbg:suggestedValue": {class: File,
      path: 663d2bcc27374715fccd8c65, name: kf_isoform_override.tsv}}
outputs:
  annotated_protected: {type: 'File[]', outputSource: rename_protected/renamed_files}
  annotated_public: {type: 'File[]', outputSource: rename_public/renamed_files}
steps:
  prefilter_vcf:
    when: $(inputs.include_expression != null)
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: input_vcf
      include_expression: bcftools_prefilter_csv
      output_basename: output_basename
    out: [filtered_vcf]
  bcftools_recontig_vcf:
    when: $(inputs.chr_rename_tsv != null)
    run: ../tools/bcftools_annotate_rename_chr.cwl
    in:
      input_vcf:
        source: [prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      chr_rename_tsv: bcftools_recontig_tsv
      output_basename: output_basename
      tool_name: tool_name
    out: [bcftools_recontig_vcf]
  normalize_vcf:
    when: $(inputs.disable_norm == false)
    run: ../tools/normalize_vcf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf:
        source: [bcftools_recontig_vcf/bcftools_recontig_vcf, prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
      disable_norm: disable_norm
    out: [normalized_vcf]
  bcftools_strip_info:
    when: $(inputs.strip_info != null)
    run: ../tools/bcftools_strip_ann.cwl
    in:
      input_vcf:
        source: [normalize_vcf/normalized_vcf, bcftools_recontig_vcf/bcftools_recontig_vcf, prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
      strip_info: bcftools_strip_columns
    out: [stripped_vcf]
  add_standard_fields:
    run: ../tools/add_strelka2_fields.cwl
    when: $(inputs.run_tool_flag)
    in:
      strelka2_vcf:
        source: [bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf, bcftools_recontig_vcf/bcftools_recontig_vcf, prefilter_vcf/filtered_vcf,
          input_vcf]
        pickValue: first_non_null
      run_tool_flag: add_common_fields
      tumor_name: input_tumor_name
      normal_name: input_normal_name
      output_basename: output_basename
    out: [output]
  vep_annotate_vcf:
    when: $(inputs.disable_annotation == false)
    run: ../tools/variant_effect_predictor_105.cwl
    in:
      reference: indexed_reference_fasta
      disable_annotation: disable_vep_annotation
      cores: vep_cores
      ram: vep_ram
      buffer_size: vep_buffer_size
      input_vcf:
        source: [add_standard_fields/output, bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf, bcftools_recontig_vcf/bcftools_recontig_vcf,
          prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
      cache: vep_cache
      merged: merged
      run_cache_existing: run_cache_existing
      run_cache_af: run_cache_af
      run_stats: run_stats
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
    out: [output_vcf]
  echtvar_anno_gnomad:
    when: $(inputs.echtvar_zips != null)
    run: ../tools/echtvar_anno.cwl
    in:
      input_vcf:
        source: [vep_annotate_vcf/output_vcf, add_standard_fields/output, bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf,
          bcftools_recontig_vcf/bcftools_recontig_vcf, prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      echtvar_zips: echtvar_anno_zips
      tbi:
        valueFrom: |
          $(1 == 1)
      output_filename:
        source: [output_basename, tool_name]
        valueFrom: |
          $(self[0]).$(self[1]).echtvar_annotated.vcf.gz
    out: [annotated_vcf]
  gatk_add_soft_filter:
    run: ../tools/gatk_variant_filter.cwl
    in:
      input_vcf:
        source: [echtvar_anno_gnomad/annotated_vcf, vep_annotate_vcf/output_vcf, add_standard_fields/output, bcftools_strip_info/stripped_vcf,
          normalize_vcf/normalized_vcf, bcftools_recontig_vcf/bcftools_recontig_vcf, prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      reference: indexed_reference_fasta
      filter_name: gatk_filter_name
      filter_expression: gatk_filter_expression
      output_basename: output_basename
      tool_name: tool_name
    out: [gatk_soft_filtered_vcf]
  hotspots_annotation:
    run: ../tools/hotspots_annotation.cwl
    in:
      input_vcf: gatk_add_soft_filter/gatk_soft_filtered_vcf
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snvs: protein_snv_hotspots
      protein_indels: protein_indel_hotspots
      output_basename: output_basename
    out: [hotspots_vcf]
  kfdrc_vcf2maf_protected:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: indexed_reference_fasta
      input_vcf: hotspots_annotation/hotspots_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      maf_center: maf_center
      custom_enst: custom_enst
    out: [output_maf]
  hard_filter_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: hotspots_annotation/hotspots_vcf
      include_expression: bcftools_public_filter
      output_basename: output_basename
      output_type:
        valueFrom: "z"
    out: [filtered_vcf]
  kfdrc_vcf2maf_public:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: indexed_reference_fasta
      input_vcf: hard_filter_vcf/filtered_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      maf_center: maf_center
      custom_enst: custom_enst
    out: [output_maf]
  rename_protected:
    run: ../tools/generic_rename_outputs.cwl
    label: Rename Protected Outputs
    in:
      input_files:
        source: [hotspots_annotation/hotspots_vcf, kfdrc_vcf2maf_protected/output_maf]
        valueFrom: |
          $([self[0], self[0].secondaryFiles[0], self[1]])
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: |
          ${
            var prefix = self[0] + '.' + self[1] + '.norm.annot.protected';
            return [
              prefix + '.vcf.gz',
              prefix + '.vcf.gz.tbi',
              prefix + '.maf'
            ];
          }
    out: [renamed_files]
  rename_public:
    run: ../tools/generic_rename_outputs.cwl
    label: Rename Public Outputs
    in:
      input_files:
        source: [hard_filter_vcf/filtered_vcf, kfdrc_vcf2maf_public/output_maf]
        valueFrom: |
          $([self[0], self[0].secondaryFiles[0], self[1]])
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: |
          ${
            var prefix = self[0] + '.' + self[1] + '.norm.annot.public';
            return [
              prefix + '.vcf.gz',
              prefix + '.vcf.gz.tbi',
              prefix + '.maf'
            ];
          }
    out: [renamed_files]
$namespaces:
  sbg: https://sevenbridges.com
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:links":
- id: 'https://github.com/kids-first/kf-annotation-tools/releases/tag/v1.3.0'
  label: github-release
