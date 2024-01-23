# Custom gnomAD Reference Creation
This documentation outlines how we created our gomAD annotation reference.
While the gnomAD source itself is very large and comprehensive, this reference is meant to add the "bare minimum" as a useful start for variant filtering out-of-the-box with a called VCF.
The main steps for reference creation are:
1. Download and pipe to bcftools and vt to grab desired fields and to normalize calls (see https://genome.sph.umich.edu/wiki/Vt#Normalization)
1. Use pysam to add custom fields including an INFO field for the `FILTER` value from the source file, a popmax AF for non-cancer populations with no bottlenecks (as defined [here](https://gnomad.broadinstitute.org/help/faf)) and a popmax __with__ bottleneck populations
1. Create an echtvar reference for blazing fast annotation of VCF files with this custom reference

## Download and normalize
Desired fields to subset from gnomAD are as follows:
```
AC
AN
AF
nhomalt
AC_popmax
AN_popmax
AF_popmax
nhomalt_popmax
AC_controls_and_biobanks
AN_controls_and_biobanks
AF_controls_and_biobanks
AF_non_cancer
primate_ai_score
splice_ai_consequence
AF_non_cancer_afr
AF_non_cancer_ami
AF_non_cancer_asj
AF_non_cancer_eas
AF_non_cancer_fin
AF_non_cancer_mid
AF_non_cancer_nfe
AF_non_cancer_oth
AF_non_cancer_raw
AF_non_cancer_sas
AF_non_cancer_amr
```
Example command:
```sh
curl -sL https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr16.vcf.bgz | bcftools annotate --threads 4  -x '^INFO/AF_non_cancer,^INFO/AF_non_cancer_afr,^INFO/AF_non_cancer_ami,^INFO/AF_non_cancer_asj,^INFO/AF_non_cancer_eas,^INFO/AF_non_cancer_fin,^INFO/AF_non_cancer_mid,^INFO/AF_non_cancer_nfe,^INFO/AF_non_cancer_oth,^INFO/AF_non_cancer_raw,^INFO/AF_non_cancer_sas,^INFO/AF_non_cancer_amr,^INFO/AC,^INFO/AN,^INFO/AF,^INFO/nhomalt,^INFO/AC_popmax,^INFO/AN_popmax,^INFO/AF_popmax,^INFO/nhomalt_popmax,^INFO/AC_controls_and_biobanks,^INFO/AN_controls_and_biobanks,^INFO/AF_controls_and_biobanks,^INFO/AF_non_cancer,^INFO/primate_ai_score,^INFO/splice_ai_consequence' | /vt/vt normalize - -n -r Homo_sapiens_assembly38.fasta  | bgzip -@4 -c > gnomad.genomes.v3.1.1.sites.chr16.bcftools_INFO_subset.vt_norm.vcf.gz
```
## Add `GNOMAD_FILTER`, `AF_non_cancer_popmax`, and `AF_non_cancer_all_popmax` INFO fields
Use `scripts/custom_vcf_info.py` to add new custom fields, calculated as follows:
 - `GNOMAD_FILTER` = `FILTER` value
 - `AF_non_cancer_popmax` max value of:
   ```
   AF_non_cancer_afr
   AF_non_cancer_amr
   AF_non_cancer_eas
   AF_non_cancer_nfe
   AF_non_cancer_sas
   ```
   when available
 - `AF_non_cancer_all_popmax` get max value of bottleneck populations:
   ```
   AF_non_cancer_ami
   AF_non_cancer_asj
   AF_non_cancer_fin
   AF_non_cancer_mid
   AF_non_cancer_oth
   ```
   and set to the greater of calculated `AF_non_cancer_popmax` or max of bottleneck populations.
Example command:
```sh
python3 custom_vcf_info.py --input_vcf gnomad.genomes.v3.1.1.sites.chr16.bcftools_INFO_subset.vt_norm.vcf.gz --output_basename gnomad.genomes.v3.1.1.sites.chr16.custom --threads 2
```
## Create echtvar reference
First, create a config json file. See [here](https://github.com/brentp/echtvar#configuration-file-for-encode) for generic details. See [docs/gnomad_update.json] for config used. This config will preppend `gnomad_3_1_1_` to all field names for source clarity upon annotation.
```sh
echtvar \
   encode \
   gnomad.v3.1.1.chr16_custom.echtvar.zip \
   gnomad_update.json \
   gnomad.genomes.v3.1.1.sites.chr16.custom.INFO_added.vcf.gz
```