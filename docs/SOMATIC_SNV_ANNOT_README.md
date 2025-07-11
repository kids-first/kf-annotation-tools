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
