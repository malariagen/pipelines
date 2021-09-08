import sys
import allel
import zarr
import zipfile
import os
import shutil


def zip_zarr(zarr_path, del_orig=False):
    """
    Assembles the entire Zarr directory into an uncompressed zip file name {zarr_path}.zip.

    Parameters
    ----------
    zarr_path : str
        Path to zarr folder.
    del_orig : bool
        Whether the original zarr directory should be deleted.

    """

    output_zip_path = zarr_path + ".zip"
    
    with zipfile.ZipFile(file=output_zip_path, mode="w", compression=zipfile.ZIP_STORED,
                         allowZip64=False) as zh:
        for curr_dir_path, dir_basenames, file_basenames in os.walk(top=zarr_path, topdown=True):
            for file_basename in file_basenames:
                # We need the absolute path to add files into the zipped zarr
                # but we want the zipped files named using relative paths to the zarr directory.
                curr_dir_rel_path = os.path.relpath(curr_dir_path, start=zarr_path)
                file_rel_path = os.path.join(curr_dir_rel_path, file_basename)
                file_abs_path = os.path.join(curr_dir_path, file_basename)
                zh.write(filename=file_abs_path, arcname=file_rel_path)

    # reopen and test
    with zipfile.ZipFile(file=output_zip_path, mode="r") as zh:
        bad = zh.testzip()
        if bad:
            raise RuntimeError(f"zip test failed, first bad file: {bad}")

    # clean up
    if del_orig:
        shutil.rmtree(zarr_path)


def main():

    import argparse
    parser = argparse.ArgumentParser(description="Convert a single sample VCF to Zarr using "
                                                 "Blosc compression")
    parser.add_argument("--input",
                        required=True,
                        help="path to input VCF file.")
    parser.add_argument("--output",
                        required=True,
                        help="path to output Zarr directory.")
    parser.add_argument("--sample",
                        required=True,
                        help="Sample identifier.")
    parser.add_argument("--contig",
                        required=False,
                        action='append',
                        dest='contigs',
                        help="Contig to extract. Multiple values may be provided. Reads VCF contigs by default.")
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
    parser.add_argument("--zip", action="store_true",
                        help="If flag exists, entire zarr folder is zipped with no "
                             "compression, and the original zarr folder is deleted. "
                             "The zip file name is the value given for the --output "
                             "argument, appended with '.zip'.")

    args = parser.parse_args()
    input_vcf_path = args.input
    output_zarr_path = args.output
    sample = args.sample
    compress_algo = args.compress_algo
    compress_level = args.compress_level
    compress_shuffle = args.compress_shuffle
    alt_number = args.alt_number
    tabix = args.tabix
    chunk_length = args.chunk_length
    chunk_width = args.chunk_width
    do_zip = args.zip
    contigs = args.contigs or [] # If no contigs provided, read from vcf file
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
    
    if not contigs:
        if input_vcf_path.endswith((".gz", ".bgz")):
            vcf_opener = gzip.open
        else:
            vcf_opener = open

        with vcf_opener(input_vcf_path, mode="rt") as vcf_open:
            for line in vcf_open:
                # Assuming vcf header follows standard vcf convention
                if line.startswith("#"):
                    if "contig=<ID=" in line:
                        chrom = line.split("<ID=")[1].split(",")[0]
                        contigs.append(chrom)
                else:
                    # Header is finished, no need to read the rest
                    break
               
    try:
        if not contigs:
            # If there still aren't any contigs present, the vcf was faulty
            log.write(f"contigs argument is empty and provided VCF does not have a valid contig header line")
            sys.exit(1)
            
        for contig in contigs:

            allel.vcf_to_zarr(
                input=input_vcf_path,
                output=output_zarr_path,
                group=f"{sample}/{contig}",
                region=contig,
                samples=[sample],
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

    if do_zip:
        zip_zarr(zarr_path=output_zarr_path, del_orig=True)


if __name__ == "__main__":
    main()
