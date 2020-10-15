from datetime import datetime
from os.path import abspath
import argparse
import numpy as np
import zarr
import allel


def encode_genotype(g):
    return '/'.join([encode_allele(a) for a in g])


def encode_allele(a):
    if a >= 0:
        return str(a)
    else:
        return "."


def main():

    parser = argparse.ArgumentParser(
        description="Subset genotypes to biallelic sites for phasing.")
    parser.add_argument("--sample-genotypes",
                        required=True,
                        help="Path to input zipped zarr file with sample genotypes.")
    parser.add_argument("--sites-called",
                        required=True,
                        help="Path to input zipped zarr file with sites at which genotypes were "
                             "called.")
    parser.add_argument("--sites-selected",
                        required=True,
                        help="Path to input zipped zarr file with sites to subset to.")
    parser.add_argument("--output",
                        required=True,
                        help="Path to output VCF file.")
    parser.add_argument("--allow-half-missing", action="store_true")
    parser.add_argument("--contig",
                        required=True,
                        action="append",
                        dest="contigs",
                        help="Contig to extract. Multiple values may be provided.")

    args = parser.parse_args()
    sites_called_path = abspath(args.sites_called)
    sites_selected_path = abspath(args.sites_selected)
    sample_genotypes_path = abspath(args.sample_genotypes)
    output_path = abspath(args.output)
    now = datetime.now().isoformat()

    # open inputs
    sites_called = zarr.open(sites_called_path, mode='r')
    sites_selected = zarr.open(sites_selected_path, mode='r')
    sample_genotypes = zarr.open(sample_genotypes_path, mode='r')

    # discover sample identifier, assume it's the top group in the sample_genotypes zarr hierarchy
    sample_id = list(sample_genotypes)[0]

    with open(output_path, mode="w", newline="\n") as out:

        # write headers
        print("##fileformat=VCFv4.3", file=out)
        print(f"##fileDate={now}", file=out)
        print("##source=malariagen/pipelines/scripts/sample_select_variants", file=out)
        print(f"##sites_called={sites_called_path}", file=out)
        print(f"##sites_selected={sites_selected_path}", file=out)
        print(f"##sample_genotypes={sample_genotypes_path}", file=out)
        print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=out)
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}", file=out)

        for contig in args.contigs:

            # load called sites data
            source_pos = allel.SortedIndex(sites_called[contig]['variants/POS'][:])
            source_ref = sites_called[contig]['variants/REF'][:]
            source_alt = sites_called[contig]['variants/ALT'][:]

            # load selected sites data
            dest_pos = sites_selected[contig]['variants/POS'][:]
            dest_ref = sites_selected[contig]['variants/REF'][:]
            dest_alt = sites_selected[contig]['variants/ALT'][:]
            dest_alleles = np.concatenate([dest_ref[:, None], dest_alt], axis=1)

            # load genotypes
            source_gt = allel.GenotypeArray(sample_genotypes[sample_id][contig]['calldata/GT'][:])

            # select sites
            loc_subset = source_pos.locate_keys(dest_pos)
            source_gt_subset = source_gt.compress(loc_subset, axis=0)
            source_ref_subset = source_ref.compress(loc_subset, axis=0)
            source_alt_subset = source_alt.compress(loc_subset, axis=0)
            assert (
                source_gt_subset.shape[0] ==
                source_ref_subset.shape[0] ==
                source_alt_subset.shape[0] ==
                dest_alleles.shape[0]
            )

            # recode alleles
            mapping = allel.create_allele_mapping(source_ref_subset, source_alt_subset,
                                                  dest_alleles)
            dest_gt = source_gt_subset.map_alleles(mapping, copy=False)

            # deal with half-missing
            if not args.allow_half_missing:
                loc_missing = np.any(dest_gt < 0, axis=2)
                dest_gt[loc_missing] = -1

            # write out
            for p, r, a, g in zip(dest_pos, dest_ref, dest_alt, dest_gt[:, 0]):
                row = [contig,  # CHROM
                       str(p),  # POS
                       ".",  # ID
                       r,  # REF
                       a[0],  # ALT
                       ".",  # QUAL
                       "PASS",  # FILTER
                       ".",  # INFO
                       "GT",  # FORMAT
                       encode_genotype(g)]
                print("\t".join(row), file=out)


if __name__ == "__main__":

    main()
