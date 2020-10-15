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
    parser.add_argument("--input",
                        required=True,
                        help="Path to input zipped zarr file with sample genotypes.")
    parser.add_argument("--sites",
                        required=True,
                        help="Path to input zipped zarr file with sites at which genotypes were "
                             "called.")
    parser.add_argument("--subset-sites",
                        required=True,
                        help="Path to input zipped zarr file with sites to subset to.")
    parser.add_argument("--output",
                        required=True,
                        help="Path to output VCF file.")
    parser.add_argument("--sample",
                        required=True,
                        help="Sample identifier.")
    parser.add_argument("--allow-half-missing", action="store_true")
    parser.add_argument("--contig",
                        required=True,
                        action="append",
                        dest="contigs",
                        help="Contig to extract. Multiple values may be provided.")

    args = parser.parse_args()
    with open(args.output, mode="wt", newline="\n") as out:

        print("#@TODO VCF headers", file=out)
        sites = zarr.open(args.sites, mode='r')
        subset_sites = zarr.open(args.subset_sites, mode='r')
        genotypes = zarr.open(args.input, mode='r')

        for contig in args.contigs:

            # load sites data
            source_pos = allel.SortedIndex(sites[contig]['variants/POS'][:])
            source_ref = sites[contig]['variants/REF'][:]
            source_alt = sites[contig]['variants/ALT'][:]

            # load subset sites data
            dest_pos = subset_sites[contig]['variants/POS'][:]
            dest_ref = subset_sites[contig]['variants/REF'][:]
            dest_alt = subset_sites[contig]['variants/ALT'][:]
            dest_alleles = np.concatenate([dest_ref[:, None], dest_alt], axis=1)

            # load genotypes
            source_gt = allel.GenotypeArray(genotypes[args.sample][contig]['calldata/GT'][:])

            # select sites
            loc_subset = source_pos.locate_keys(dest_pos)
            source_gt_subset = source_gt.compress(loc_subset, axis=0)
            source_ref_subset = source_ref.compress(loc_subset, axis=0)
            source_alt_subset = source_alt.compress(loc_subset, axis=0)

            # recode alleles
            mapping = allel.create_allele_mapping(source_ref_subset, source_alt_subset,
                                                  dest_alleles)
            dest_gt = source_gt_subset.map_alleles(mapping)

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
