import sys
import allel
import zarr


def main():

    import argparse
    parser = argparse.ArgumentParser(description="Convert a multiple sample VCF to Zarr using "
                                                 "Blosc compression")
    parser.add_argument("--input",
                        required=True,
                        help="path to input VCF file.")
    parser.add_argument("--output",
                        required=True,
                        help="path to output Zarr directory.")
    parser.add_argument("--contig",
                        required=True,
                        action='append',
                        dest='contigs',
                        help="Contig to extract. Multiple values may be provided.")
    parser.add_argument("--field",
                        required=True,
                        action='append',
                        dest='fields',
                        help="Field to extract, e.g., 'variants/MQ' or 'calldata/GT'. Multiple "
                             "values may be provided.")
    parser.add_argument("--compress-algo",
                        help="Blosc compression algorithm.  Choose from [zstd, blosclz, lz4, "
                             "lz4hc, zlib, snappy].",
                        default="zstd")
    parser.add_argument("--compress-level",
                        help="Compression level. Choose integer from [0, 9].",
                        default=1,
                        type=int)
    parser.add_argument("--compress-shuffle",
                        help="Type of data shuffling used to obtain contiguous runs of same "
                             "values for improving compression. Choose integer value from "
                             "NOSHUFFLE (0), SHUFFLE (1), BITSHUFFLE (2) or AUTOSHUFFLE (-1). "
                             "If -1 (default), bit-shuffle will be used for buffers with "
                             "itemsize 1, and byte-shuffle will be used otherwise.",
                        default=-1,
                        type=int)
    parser.add_argument("--dtype",
                        help="TODO",
                        action="append",
                        type=str)
    parser.add_argument("--alt-number",
                        help="Expected maximum number of alternate alleles.",
                        default=3,
                        type=int)
    parser.add_argument("--tabix",
                        help="Path to tabix executable v0.2.5+.",
                        default="tabix")
    parser.add_argument("--chunk-length",
                        help="Chunk length in number of variants.",
                        default=2**18,
                        type=int)
    parser.add_argument("--chunk-width",
                        help="Chunk width in number of samples.",
                        default=64,
                        type=int)
    parser.add_argument("--log",
                        help="Path to logfile, stdout or stderr. Default: stdout.",
                        default="stdout")

    args = parser.parse_args()
    input_vcf_path = args.input
    output_zarr_path = args.output
    compress_algo = args.compress_algo
    compress_level = args.compress_level
    compress_shuffle = args.compress_shuffle
    alt_number = args.alt_number
    tabix = args.tabix
    chunk_length = args.chunk_length
    chunk_width = args.chunk_width
    contigs = args.contigs
    fields = args.fields
    log = args.log.strip()

    log_file_needs_closing = False
    if log == "stderr":
        log_file = sys.stderr
    elif log == "stdout":
        log_file = sys.stdout
    else:
        log_file = open(log, "w")
        log_file_needs_closing = True

    # TODO add dtype handling

    try:

        for contig in contigs:

            allel.vcf_to_zarr(
                input=input_vcf_path,
                output=output_zarr_path,
                group=contig,
                region=contig,
                compressor=zarr.Blosc(cname=compress_algo, clevel=compress_level,
                                      shuffle=compress_shuffle),
                overwrite=True,
                tabix=tabix,
                fields=fields,
                alt_number=alt_number,
                chunk_length=chunk_length,
                chunk_width=chunk_width,
                log=log_file,
            )

    finally:
        if log_file_needs_closing:
            log_file.close()


if __name__ == "__main__":
    main()
