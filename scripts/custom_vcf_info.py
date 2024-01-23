"""
Add custom INFO fields to vcf using pysam
"""
import pysam
import argparse


def create_mod_vcf(output_path, input_path, threads):
    """ Create a new VCF in which standard FORMAT fields have been
            calculated and added=based on those provided natively by Strelka2

        Args:
            output_path (str): path to the VCF being created
            input_path (str): path to the Strelka2 VCF
            tumor_name (str): name of tumor sample in input
            normal_name (str) name of normal sample in input

        Raises:
            IOError if a sample name in the input is neither tumor_name
                nor normal_name
    """

    input_vcf = pysam.VariantFile(input_path, 'r', threads=threads)

    input_vcf.header.info.add('GNOMAD_FILTER', '.', 'String',
            'Value of FILTER for gnomAD variant. Use to include/exclude non-PASS variants')
    input_vcf.header.info.add('AF_non_cancer_popmax', 'A', 'Float',
            ('Max AF of non-bottleneck populations in AF_non_cancer'))
    input_vcf.header.info.add('AF_non_cancer_all_popmax', 'A', 'Float',
            ('Max AF of populations in AF_non_cancer INCLUDING bottleneck'))

    output = pysam.VariantFile(output_path, 'w', header=input_vcf.header)
    # non-bottleneck population fields to use for popmax
    pop_fields = ['AF_non_cancer_afr', 'AF_non_cancer_amr', 'AF_non_cancer_eas', 'AF_non_cancer_nfe', 'AF_non_cancer_sas']
    # bottleneck population fields
    pop_fields_bn = ['AF_non_cancer_ami', 'AF_non_cancer_asj', 'AF_non_cancer_fin', 'AF_non_cancer_mid', 'AF_non_cancer_oth']
    for record in input_vcf.fetch():
        try:
            record.info['GNOMAD_FILTER'] = record.filter.keys()
            if 'AF_non_cancer' in record.info:
                pop_values = []
                pop_values_bn = []
                pop_values = [record.info[pop] for pop in pop_fields if pop in record.info]
                # get non bottleneck pop max
                if len(pop_values) > 0:
                    record.info['AF_non_cancer_popmax'] = max(pop_values)
                # get pop max with bottleneck, just compare to non
                pop_values_bn = [record.info[pop] for pop in pop_fields_bn if pop in record.info]
                if len(pop_values_bn) > 0:
                    if len(pop_values) > 0:
                        record.info['AF_non_cancer_all_popmax'] = max(max(pop_values), max(pop_values_bn))
                    else:
                        record.info['AF_non_cancer_all_popmax'] = max(pop_values_bn)
            output.write(record)
        except Exception as e:
            print(e)
            print("Failed at {} {}".format(record.contig, record.pos))
            exit(1)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
            description = 'Add custom fields to gnomAD vcf. Limited scope')

    parser.add_argument('--input_vcf', 
            help='VCF to add INFO to')
    parser.add_argument('--output_basename',
            help='String to use as basename for output file [e.g.] task ID')
    parser.add_argument('--threads',
            help='Num threads to use for read/write', default=1)

    args = parser.parse_args()

    output_vcf_name = args.output_basename + ".INFO_added.vcf.gz"

    # Create and index the modified VCF
    create_mod_vcf(output_vcf_name, args.input_vcf, int(args.threads))
    pysam.tabix_index(output_vcf_name, preset="vcf", force=True)
