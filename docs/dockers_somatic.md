# Dockers of kfdrc-somatic-snv-annot-workflow.cwl

TOOL|DOCKER
-|-
add_strelka2_fields.cwl|pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0
bcftools_annotate.cwl|pgc-images.sbgenomics.com/d3b-bixu/bcftools:1.20
bcftools_filter_vcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/bcftools:1.20
bcftools_norm.cwl|pgc-images.sbgenomics.com/d3b-bixu/bcftools:1.20
echtvar_anno.cwl|pgc-images.sbgenomics.com/d3b-bixu/echtvar:0.2.0
gatk_variant_filter.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
generic_rename_outputs.cwl|None
hotspots_annotation.cwl|quay.io/biocontainers/pysam:0.21.0--py310h41dec4a_1
kf_mskcc_vcf2maf.cwl|pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3
variant_effect_predictor_105.cwl|ensemblorg/ensembl-vep:release_105.0
