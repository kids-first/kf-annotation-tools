cwlVersion: v1.2
class: Workflow
id: kfdrc-germline-snv-annot-wf
label: Kids First DRC Germline SNV Annotation Workflow

doc: |-
  # Kids First DRC Germline SNV Annotation Workflow
  This workflow is used to annotate germline outputs with popular annotation resources. This includes using VEP to annotate with ENSEMBL v105 reference as well using bcftools to add further annotation described below.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  ## Overall annotation steps
  1. Prefilter input VCF (optional) to remove variants that are undesired to go into annotation
  1. Normalize VCF
  1. Strip pre-existing annotations (optional) to prevent downstream conflicts
  1. Annotate with VEP 105. Optional plugins include:
     - dbnsfp
     - cadd
  1. Use echtvar to annotate with an external reference (default gnomad 3.1.1)
  1. Use bcftools to annotate with another external reference (optional clinvar)
  1. Simple rename outputs step

  ## Default annotations
  By default, the workflow will add the following annotations:

  ### ENSEMBL 105
  This is added on using variant effect predictor to use the ENSEMBL reference to add gene model information as well as additional resources provided in their cache. It's highly recommended that when you download their cache, to [convert and index](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#convert). It will speed up annotation and reduce memory footprint significantly. Annotation resources in the cache include:

  ```
  # CACHE UPDATED 2022-09-26 18:18:29
  assembly	GRCh38
  bam	GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bam
  polyphen	b
  sift	b
  source_assembly	GRCh38.p13
  source_gencode	GENCODE 39
  source_genebuild	2014-07
  source_polyphen	2.2.2
  source_refseq	2021-05-28 21:42:08 - GCF_000001405.39_GRCh38.p13_genomic.gff
  source_sift	sift5.2.2
  species	homo_sapiens
  variation_cols	chr,variation_name,failed,somatic,start,end,allele_string,strand,minor_allele,minor_allele_freq,clin_sig,phenotype_or_disease,clin_sig_allele,pubmed,var_synonyms,AFR,AMR,EAS,EUR,SAS,AA,EA,gnomAD,gnomAD_AFR,gnomAD_AMR,gnomAD_ASJ,gnomAD_EAS,gnomAD_FIN,gnomAD_NFE,gnomAD_OTH,gnomAD_SAS
  source_COSMIC	94
  source_HGMD-PUBLIC	20204
  source_ClinVar	105202106
  source_dbSNP	154
  source_1000genomes	phase3
  source_ESP	V2-SSA137
  source_gnomAD	r2.1.1
  regulatory	1
  cell_types	A549,A673,B,B_(PB),CD14+_monocyte_(PB),CD14+_monocyte_1,CD4+_CD25+_ab_Treg_(PB),CD4+_ab_T,CD4+_ab_T_(PB)_1,CD4+_ab_T_(PB)_2,CD4+_ab_T_(Th),CD4+_ab_T_(VB),CD8+_ab_T_(CB),CD8+_ab_T_(PB),CMP_CD4+_1,CMP_CD4+_2,CMP_CD4+_3,CM_CD4+_ab_T_(VB),DND-41,EB_(CB),EM_CD4+_ab_T_(PB),EM_CD8+_ab_T_(VB),EPC_(VB),GM12878,H1-hESC_2,H1-hESC_3,H9_1,HCT116,HSMM,HUES48,HUES6,HUES64,HUVEC,HUVEC-prol_(CB),HeLa-S3,HepG2,K562,M0_(CB),M0_(VB),M1_(CB),M1_(VB),M2_(CB),M2_(VB),MCF-7,MM.1S,MSC,MSC_(VB),NHLF,NK_(PB),NPC_1,NPC_2,NPC_3,PC-3,PC-9,SK-N.,T_(PB),Th17,UCSF-4,adrenal_gland,aorta,astrocyte,bipolar_neuron,brain_1,cardiac_muscle,dermal_fibroblast,endodermal,eosinophil_(VB),esophagus,foreskin_fibroblast_2,foreskin_keratinocyte_1,foreskin_keratinocyte_2,foreskin_melanocyte_1,foreskin_melanocyte_2,germinal_matrix,heart,hepatocyte,iPS-15b,iPS-20b,iPS_DF_19.11,iPS_DF_6.9,keratinocyte,kidney,large_intestine,left_ventricle,leg_muscle,lung_1,lung_2,mammary_epithelial_1,mammary_epithelial_2,mammary_myoepithelial,monocyte_(CB),monocyte_(VB),mononuclear_(PB),myotube,naive_B_(VB),neuron,neurosphere_(C),neurosphere_(GE),neutro_myelocyte,neutrophil_(CB),neutrophil_(VB),osteoblast,ovary,pancreas,placenta,psoas_muscle,right_atrium,right_ventricle,sigmoid_colon,small_intestine_1,small_intestine_2,spleen,stomach_1,stomach_2,thymus_1,thymus_2,trophoblast,trunk_muscle
  source_regbuild	1.0
  var_type	tabix
  ```

  ### [gnomAD 3.1.1](https://gnomad.broadinstitute.org/)
  Using echtvar, we annotate from a [custom implementation of gnomAD v3.1.1](CUSTOM_GNOMAD_REF.md) the following population statistics (columns are give a `gnomad_3_1_1_` prefix to denote source):
  ```
  gnomad_3_1_1_AC
  gnomad_3_1_1_AN
  gnomad_3_1_1_AF
  gnomad_3_1_1_nhomalt
  gnomad_3_1_1_AC_popmax
  gnomad_3_1_1_AN_popmax
  gnomad_3_1_1_AF_popmax
  gnomad_3_1_1_nhomalt_popmax
  gnomad_3_1_1_AC_controls_and_biobanks
  gnomad_3_1_1_AN_controls_and_biobanks
  gnomad_3_1_1_AF_controls_and_biobanks
  gnomad_3_1_1_AF_non_cancer
  gnomad_3_1_1_primate_ai_score
  gnomad_3_1_1_splice_ai_consequence
  gnomad_3_1_1_AF_non_cancer_afr
  gnomad_3_1_1_AF_non_cancer_ami
  gnomad_3_1_1_AF_non_cancer_asj
  gnomad_3_1_1_AF_non_cancer_eas
  gnomad_3_1_1_AF_non_cancer_fin
  gnomad_3_1_1_AF_non_cancer_mid
  gnomad_3_1_1_AF_non_cancer_nfe
  gnomad_3_1_1_AF_non_cancer_oth
  gnomad_3_1_1_AF_non_cancer_raw
  gnomad_3_1_1_AF_non_cancer_sas
  gnomad_3_1_1_AF_non_cancer_amr
  gnomad_3_1_1_AF_non_cancer_popmax
  gnomad_3_1_1_AF_non_cancer_all_popmax
  gnomad_3_1_1_FILTER
  ```

  ## Optional annotations
  ### [dbNSFP v4.3a](http://database.liulab.science/dbNSFP#intro)
  This resource compiles from dozens of sources annotations for ~84M SNVs. By default, from this resource, we recommend the following:
  ```
  SIFT4G_pred
  Polyphen2_HDIV_pred
  Polyphen2_HVAR_pred
  LRT_pred
  MutationTaster_pred
  MutationAssessor_pred
  FATHMM_pred
  PROVEAN_pred
  VEST4_score
  VEST4_rankscore
  MetaSVM_pred
  MetaLR_pred
  MetaRNN_pred
  M-CAP_pred
  REVEL_score
  REVEL_rankscore
  PrimateAI_pred
  DEOGEN2_pred
  BayesDel_noAF_pred
  ClinPred_pred
  LIST-S2_pred
  Aloft_pred
  fathmm-MKL_coding_pred
  fathmm-XF_coding_pred
  Eigen-phred_coding
  Eigen-PC-phred_coding
  phyloP100way_vertebrate
  phyloP100way_vertebrate_rankscore
  phastCons100way_vertebrate
  phastCons100way_vertebrate_rankscore
  TWINSUK_AC
  TWINSUK_AF
  ALSPAC_AC
  ALSPAC_AF
  UK10K_AC
  UK10K_AF
  gnomAD_exomes_controls_AC
  gnomAD_exomes_controls_AN
  gnomAD_exomes_controls_AF
  gnomAD_exomes_controls_nhomalt
  gnomAD_exomes_controls_POPMAX_AC
  gnomAD_exomes_controls_POPMAX_AN
  gnomAD_exomes_controls_POPMAX_AF
  gnomAD_exomes_controls_POPMAX_nhomalt
  Interpro_domain
  GTEx_V8_gene
  GTEx_V8_tissue
  ```

  ### [CADD v1.6](https://cadd.gs.washington.edu/)
  Using a VEP plugin, we add Combined Annotation Dependent Depletion scores

  ### [ClinVar 20220507](https://www.ncbi.nlm.nih.gov/clinvar/)
  A curated resource with annotations of clinical significance per variant. Note, for this pipeline, the default reference was modified by:
     - Switching from `1` chromosome nomenclature to `chr1`, and especially `MT` -> `chrM`
     - Removing the entry assigned to `NW_009646201.1`. It's a benign it and also not present in our fasta reference.
  We recommend the following:
  ```
  ALLELEID
  CLNDN
  CLNDNINCL
  CLNDISDB
  CLNDISDBINCL
  CLNHGVS
  CLNREVSTAT
  CLNSIG
  CLNSIGCONF
  CLNSIGINCL
  CLNVC
  CLNVCSO
  CLNVI
  ```

  ### [InterVar](https://github.com/WGLab/InterVar)
  This is a custom reference generated by the authors of the tool linked above. It contains only exonic snps. To utilize the full capabilities of their classification, you must run the tool.

  ## Workflow Inputs
  ### Required
   - `indexed_reference_fasta` files: Homo_sapiens_assembly38.fasta, Homo_sapiens_assembly38.dict, Homo_sapiens_assembly38.fasta.fai
   - `input_vcf`: Input vcf file to annotate
   - `output_basename`: string prefix of outputs
   - `tool_name`: short descriptive string of tool output being annotated
  ### RECOMMENDED
   - `echtvar_anno_zips` file array: Annotation ZIP files for echtvar anno
   - `vep_cache` file: `homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz`
   - `merged` boolean: Set to true if merged cache used, default: `true`
   - `run_cache_existing` boolean: Run the check_existing flag for cache, default: `true`
   - `run_cache_af` boolean: Run the allele frequency flags for cache, default: `true`
   - `run_stats` boolean: Create stats file. Disable for speed, default: `false`

  ### Optional
   - `bcftools_prefilter_csv`: CSV of bcftools filter params if you want to prefilter before annotation
   - `disable_normalization` boolean: Skip normalizing if input is already normed, default is `false`
   - `bcftools_strip_columns`: CSV string of columns to strip if needed to avoid conflict, i.e INFO/AF
   - `vep_ram` int: In GB, may need to increase this value depending on the size/complexity of input, default: `48`
   - `vep_cores` int: Number of cores to use. May need to increase for really large inputs, default: `32`,
   - `vep_buffer_size` int: Increase or decrease to balance speed and memory usage, default: `100000`
   - `bcftools_annot_clinvar_columns`: CSV string of columns from annotation to port into the input VCF, default: `INFO/ALLELEID,INFO/CLNDN,INFO/CLNDNINCL,INFO/CLNDISDB,INFO/CLNDISDBINCL,INFO/CLNHGVS,INFO/CLNREVSTAT,INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNSIGINCL,INFO/CLNVC,INFO/CLNVCSO,INFO/CLNVI`
   - `clinvar_annotation_vcf` files: clinvar_20220507_chr_fixed.vcf.gz, clinvar_20220507_chr_fixed.vcf.gz.tbi
   - `dbnsfp` file: dbNSFP4.3a_grch38.gz, dbNSFP4.3a_grch38.gz.tbi, dbNSFP4.3a_grch38.readme.txt
   - `dbnsfp_fields` string: CSV string with desired fields to annotate. Use ALL to grab all, default: `SIFT4G_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_pred,MetaLR_pred,MetaRNN_pred,M-CAP_pred,REVEL_score,REVEL_rankscore,PrimateAI_pred,DEOGEN2_pred,BayesDel_noAF_pred,ClinPred_pred,LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Eigen-phred_coding,Eigen-PC-phred_coding,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,UK10K_AC,UK10K_AF,gnomAD_exomes_controls_AC,gnomAD_exomes_controls_AN,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_nhomalt,gnomAD_exomes_controls_POPMAX_AC,gnomAD_exomes_controls_POPMAX_AN,gnomAD_exomes_controls_POPMAX_AF,gnomAD_exomes_controls_POPMAX_nhomalt,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue`
   - `cadd_indels` file: CADDv1.6-38-gnomad.genomes.r3.0.indel.tsv.gz, CADDv1.6-38-gnomad.genomes.r3.0.indel.tsv.gz.tbi
   - `cadd_snvs` file: CADDv1.6-38-whole_genome_SNVs.tsv.gz, CADDv1.6-38-whole_genome_SNVs.tsv.gz.tbi
   - `intervar` file: Exons.all.hg38.intervar.2021-07-31.vcf.gz, Exons.all.hg38.intervar.2021-07-31.vcf.gz.tbi

  ## Workflow Outputs
   - `annotated_vcf` file: VCF file with all applied annotations

  ## Notes On Reference Generation
  Currently, VEP not only provides gene model annotation, but also allows for additional annotations to be added. Therefore, for CADD and dbNSFP, existing files were formatted in order to use VEP plugins. Please see [their documentation](https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html) for information on how the refernces were generated.
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict], "sbg:suggestedValue": {class: File, path: 60639014357c3a53540ca7a3,
      name: Homo_sapiens_assembly38.fasta, secondaryFiles: [{class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict},
        {class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai}]}}
  input_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Input vcf to annotate"}
  output_basename: string
  tool_name: {type: string, doc: "File name string suffx to use for output files"}

  bcftools_prefilter_csv: {type: 'string?', doc: "csv of bcftools filter params if you want to prefilter before annotation"}
  disable_normalization: {type: 'boolean?', doc: "Skip normalizing if input is already normed", default: false}
  # bcftools strip, if needed
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  # bcftools annotate if more to do
  bcftools_annot_clinvar_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf", default: "INFO/ALLELEID,INFO/CLNDN,INFO/CLNDNINCL,INFO/CLNDISDB,INFO/CLNDISDBINCL,INFO/CLNHGVS,INFO/CLNREVSTAT,INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNSIGINCL,INFO/CLNVC,INFO/CLNVCSO,INFO/CLNVI"}
  clinvar_annotation_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "additional bgzipped annotation vcf file"}
  echtvar_anno_zips: { type: 'File[]?', doc: "Annotation ZIP files for echtvar anno",
    "sbg:suggestedValue": [{class: File, path: 65c64d847dab7758206248c6, name: gnomad.v3.1.1.custom.echtvar.zip}] } 
  # VEP-specific
  disable_vep_annotation: {type: 'boolean?', doc: "Disable VEP Annotation and skip this task.", default: false}
  vep_ram: {type: 'int?', default: 48, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', default: 32, doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', default: 100000, doc: "Increase or decrease to balance speed and memory usage"}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache", "sbg:suggestedValue": {class: File, path: 6332f8e47535110eb79c794f,
      name: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz}}
  dbnsfp: {type: 'File?', secondaryFiles: [.tbi, ^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing
      dbNSFP annotations"}
  dbnsfp_fields: {type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all", default: 'SIFT4G_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_pred,MetaLR_pred,MetaRNN_pred,M-CAP_pred,REVEL_score,REVEL_rankscore,PrimateAI_pred,DEOGEN2_pred,BayesDel_noAF_pred,ClinPred_pred,LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Eigen-phred_coding,Eigen-PC-phred_coding,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,UK10K_AC,UK10K_AF,gnomAD_exomes_controls_AC,gnomAD_exomes_controls_AN,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_nhomalt,gnomAD_exomes_controls_POPMAX_AC,gnomAD_exomes_controls_POPMAX_AN,gnomAD_exomes_controls_POPMAX_AF,gnomAD_exomes_controls_POPMAX_nhomalt,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue'}
  merged: {type: 'boolean?', doc: "Set to true if merged cache used", default: true}
  run_cache_existing: {type: 'boolean?', doc: "Run the check_existing flag for cache", default: true}
  run_cache_af: {type: 'boolean?', doc: "Run the allele frequency flags for cache", default: true}
  run_stats: {type: 'boolean?', doc: "Create stats file? Disable for speed", default: false}
  cadd_indels: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations"}
  cadd_snvs: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations"}
  intervar: {type: 'File?', doc: "Intervar vcf-formatted file. Exonic SNVs only - for more comprehensive run InterVar. See docs for
      custom build instructions", secondaryFiles: [.tbi]}

outputs:
  annotated_vcf: {type: 'File[]', outputSource: rename_output/renamed_files}

steps:
  prefilter_vcf:
    when: $(inputs.include_expression != null)
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: input_vcf
      include_expression: bcftools_prefilter_csv
      output_basename: output_basename
    out: [filtered_vcf]

  normalize_vcf:
    when: $(inputs.disable_norm == false)
    run: ../tools/normalize_vcf.cwl
    in:
      disable_norm: disable_normalization
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf:
        source: [prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
    out: [normalized_vcf]

  bcftools_strip_info:
    when: $(inputs.strip_info != null)
    run: ../tools/bcftools_strip_ann.cwl
    in:
      input_vcf:
        source: [normalize_vcf/normalized_vcf, prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
      strip_info: bcftools_strip_columns
    out: [stripped_vcf]

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
        source: [bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf, prefilter_vcf/filtered_vcf, input_vcf]
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
      intervar: intervar
    out: [output_vcf]

  echtvar_anno:
    when: $(inputs.echtvar_zips != null)
    run: ../tools/echtvar_anno.cwl
    in:
      input_vcf:
        source: [vep_annotate_vcf/output_vcf, bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf, prefilter_vcf/filtered_vcf,
          input_vcf]
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

  bcftools_clinvar_annotate:
    when: $(inputs.annotation_vcf != null)
    run: ../tools/bcftools_annotate.cwl
    in:
      input_vcf:
        source: [echtvar_anno/annotated_vcf, vep_annotate_vcf/output_vcf, bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf,
          prefilter_vcf/filtered_vcf, input_vcf]
        pickValue: first_non_null
      annotation_vcf: clinvar_annotation_vcf
      columns: bcftools_annot_clinvar_columns
      output_basename: output_basename
      tool_name: tool_name
    out: [bcftools_annotated_vcf]

  rename_output:
    run: ../tools/generic_rename_outputs.cwl
    label: Rename Outputs
    in:
      input_files:
        source: [bcftools_clinvar_annotate/bcftools_annotated_vcf, echtvar_anno/annotated_vcf, vep_annotate_vcf/output_vcf]
        valueFrom: |
          ${ var first_non_null = self.filter(function(e) { return e != null }).shift();
            return [first_non_null, first_non_null.secondaryFiles[0]];
          }
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: |
          ${ var pro_vcf = [self[0], self[1], 'vcf.gz'].join('.');
            return [pro_vcf, pro_vcf + '.tbi'];
          }
    out: [renamed_files]

$namespaces:
  sbg: https://sevenbridges.com

sbg:license: Apache License 2.0
sbg:publisher: KFDRC

"sbg:links":
- id: 'https://github.com/kids-first/kids-first/kf-annotation-tools/releases/tag/v1.1.0'
  label: github-release

