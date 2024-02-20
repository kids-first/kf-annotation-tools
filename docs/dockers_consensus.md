# Dockers of kfdrc-germline-snv-annot-workflow.cwl

TOOL|DOCKER
-|-
bcftools_annotate.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
bcftools_filter_vcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest
bcftools_strip_ann.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
echtvar_anno.cwl|pgc-images.sbgenomics.com/d3b-bixu/echtvar:0.2.0
generic_rename_outputs.cwl|None
normalize_vcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
variant_effect_predictor_105.cwl|ensemblorg/ensembl-vep:release_105.0
