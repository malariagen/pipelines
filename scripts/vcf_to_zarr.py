import sys
import allel
import zarr
import zipfile
import os
import shutil



def auto_int(x):
    """
    Automatic base detection
    """
    return int(x, 0)


def is_str(s):
    """
    Returns whether s is a string.  Supports python v2 and v3.
    """
    if sys.version_info[0] == 2:
        return isinstance(s, basestring)
    elif sys.version_info[0] == 3:
        return isinstance(s, str)
    else:
        raise ValueError("Unsupported python version {}".format(sys.version_info[0]))


def get_iter(mydict):
    """
    Retrieves iterator over tuples of (key, value) in dict.  Supports python v2 and v3.
    """
    if sys.version_info[0] == 2:
        return mydict.iteritems()
    elif sys.version_info[0] == 3:
        return mydict.items()
    else:
        raise ValueError("Unsupported python version {}".format(sys.version_info[0]))


def build_sample_zarr(input_vcf_path, output_zarr_path, group, region, fields, alt_number,
                      tabix, compress_algo, compress_level, compress_shuffle,
                      chunk_length, chunk_width, log, permissions=0o660):
    """
    Wrapper method for scikit's allel.vcf_to_zarr()

    Converts a VCF to Zarr.

    If the region-filtered VCF is empty, then prints a warning to stdout and quits.

    Although tabix supports multiple regions as separate commandline parameters,
    scikit-allele (as of v1.1.9) passes in those regions as a single string (eg "region1 region2 region3") to tabix.
    In effect, only single regions are allowed, since tabix is unable to parse the concatenated region string.

    Args:
        input_vcf_path (str): path to input VCF
        group (str): Group name
        region (str or list):  str region to extract from VCF in tabix format.
                        If list of str, then each region will be extracted from
                        the VCF and placed in a subgroup under the parent group
                        specifed in group.
        fields (list):  list of fields to extract from VCF.  Must be in format used by allel
        alt_number (int):  total alternate alleles
        tabix (str):  path to tabix executable.  Only required if region specified.
        compress_algo (str):  Blosc compression algorithm.  Choose from [zstd, blosclz, lz4, lz4hc, zlib, snappy].
        compress_level (int):  Blosc compression level.  Choose integer from [0, 9].
        compress_shuffle (int):  Blosc compression shuffling type.
                        Choose integer value from NOSHUFFLE (0), SHUFFLE (1), BITSHUFFLE (2) or AUTOSHUFFLE (-1).
                        If -1 (default), bit-shuffle will be used for buffers with itemsize 1,
                        and byte-shuffle will be used otherwise.
        chunk_length (int):  chunk length in Bytes
        chunk_width (int):  chunk width in Bytes
        log (str):  path to logfile, stdout, or stderr.  Default is stdout.
        permissions (int):  octal notation of unix permissions on all files within the Zarr directory.
                        Ensure there is a leading 0o followed by 3 integer digits
                        representing user, group, and other permissions.
                        Default is 0o660: user allow read + allow write + forbid execute,
                        group allow read + allow write + forbid execute,
                        other forbid read + forbid write + forbid execute.
                        Note that directories will always have explicit allow execute for user, group, and other,
                        so that everyone can traverse the directories if they also have read permission.

    """

    # Allow read traversal of directories for everyone (i.e. user, group, other)
    PERM_ALLOW_DIR_TRAVERSE_MASK = 0o111

    output_zip_path = output_zarr_path + ".zip"

    group_to_region = {}
    if not region or is_str(region):
        group_to_region[group] = region
    else:
        for subgroup_region in region:
            subgroup = group + "/" + subgroup_region
            group_to_region[subgroup] = subgroup_region

    for subgroup, subgroup_region in get_iter(group_to_region):
        try:
            allel.vcf_to_zarr(
                input=input_vcf_path,
                output=output_zarr_path,
                group=subgroup,
                region=subgroup_region,
                compressor=zarr.Blosc(cname=compress_algo, clevel=compress_level, shuffle=compress_shuffle),
                overwrite=True,
                tabix=tabix,
                fields=fields,
                alt_number=alt_number,
                chunk_length=chunk_length,
                chunk_width=chunk_width,
                log=log
            )
        except StopIteration:
            # As of scikit-allel v1.1.9, allel.vcf_to_zarr will die with StopIteration
            # error if the vcf filtered for the desired region is empty.
            # We want to warn the user and quit nicely when that happens.
            print ("No entries in vcf {} filtered for region '{}'".format(input_vcf_path, subgroup_region))

            # allele.vcf_to_zarr() will create calldata and variants subgroups before it dies of StopIteration.
            # but create the subgroup again just in case.
            root_group = zarr.open_group(output_zarr_path, mode='a', path=subgroup)


    # As of scikit-allel v1.1.9 and zarr v2.1.4,
    # files within the Zarr directory are only read-writeable by the user since
    # they are created from tempfile.NamedTemporaryFiles.
    # See https://stackoverflow.com/questions/10541760/can-i-set-the-umask-for-tempfile-namedtemporaryfile-in-python
    # As a workaround, we recursively set explicit permissions on all files
    # within the Zarr directory.
    # We always allow read traversal on zarr subdirectories for everyone,
    # but allow the user to dictate permissions on the files within the zarr.
    dir_permissions = permissions | PERM_ALLOW_DIR_TRAVERSE_MASK
    for curr_dir_path, dir_basenames, file_basenames in os.walk(top=output_zarr_path, topdown=True):
        for dir_basename in dir_basenames:
            dir_abs_path = os.path.join(curr_dir_path, dir_basename)
            os.chmod(dir_abs_path, dir_permissions)
        for file_basename in file_basenames:
            file_abs_path = os.path.join(curr_dir_path, file_basename)
            os.chmod(file_abs_path, permissions)



def zip_zarr(zarr_path, del_orig=False):
    """
    Assembles the entire Zarr directory into an
    uncompressed zip file name {zarr_path}.zip

    Args:
        zarr_path (str):  path to zarr folder
        del_orig (bool):  whether the original zarr directory should be deleted
    """

    output_zip_path = zarr_path + ".zip"
    # Get handle for uncompressed zip file for writing
    zh = zipfile.ZipFile(file=output_zip_path, mode="w", compression=zipfile.ZIP_STORED, allowZip64=False)
    for curr_dir_path, dir_basenames, file_basenames in os.walk(top=zarr_path, topdown=True):
        for file_basename in file_basenames:
            # We need the absolute path to add files into the zipped zarr
            # but we want the zipped files named using relative paths to the zarr directory.
            curr_dir_rel_path = os.path.relpath(curr_dir_path, start=output_zarr_path)
            file_rel_path = os.path.join(curr_dir_rel_path, file_basename)
            file_abs_path = os.path.join(curr_dir_path, file_basename)
            zh.write(filename=file_abs_path, arcname=file_rel_path)

    # Check files in zip for CRC.  Unable to do this after we close.
    first_bad_file = zh.testzip()
    if first_bad_file:
        raise ValueError("CRC check failed for file {}".format(first_bad_file))

    zh.close()

    if del_orig:
        shutil.rmtree(zarr_path)


if __name__ == "__main__":

    try:
        # running under snakemake
        input_vcf = snakemake.input.vcf
        output_zarr_path = snakemake.output.zarr
        zarr_group = snakemake.params.group
        region = snakemake.params.region
        compress_algo = snakemake.params.compress_algo
        compress_level = snakemake.params.compress_level
        compress_shuffle = snakemake.params.compress_shuffle
        fields = snakemake.params.fields
        alt_number = snakemake.params.alt_number
        tabix = snakemake.params.tabix
        chunk_length = snakemake.params.chunk_length
        chunk_width = snakemake.params.chunk_width
        log = snakemake.params.log
        is_zip = snakemake.params.is_zip
        permissions = snakemake.params.permissions

    except NameError:
        # running from command line
        import argparse
        parser = argparse.ArgumentParser(description='Convert VCF to Zarr using Blosc compression')
        parser.add_argument('--input_vcf', help='path to input VCF file. Required.')
        parser.add_argument('--output_zarr_path', help='path to output Zarr directory. Required.')
        parser.add_argument('--group', help='comma separated list of Zarr hierarchy group names')
        parser.add_argument('--compress_algo',
                            help=("Blosc compression algorithm.  Choose from [zstd, blosclz, lz4, lz4hc, zlib, snappy].  " +
                                  "default: %(default)s"),
                            default="zstd")
        parser.add_argument('--compress_level',
                            help=("Compression level.  Choose integer from [0, 9].  " +
                                  "default: %(default)s"),
                            default=5,
                            type=int)
        parser.add_argument('--compress_shuffle',
                            help=("Type of data shuffling used to obtain contiguous runs of same values for improving compression.  " +
                                  "Choose integer value from NOSHUFFLE (0), SHUFFLE (1), BITSHUFFLE (2) or AUTOSHUFFLE (-1).  " +
                                  "If -1 (default), bit-shuffle will be used for buffers with itemsize 1, and byte-shuffle will be used otherwise.  " +
                                  "default: %(default)s"),
                            default=2,
                            type=int)
        parser.add_argument('--region',
                            help='Region string in tabix format. Multiple regions delimited by commas. ' +
                            'EG) chr1:2-20,chr2.  ' +
                            'If multiple regions specified, each region will be placed in a separate subgroup under --group')
        parser.add_argument('--fields', help="Comma delimited fields to extract.  default: %(default)s",
                            default="calldata/GT,calldata/GQ,calldata/DP,calldata/AD")
        parser.add_argument('--alt_number', help='Expected number of alternate alleles. default: %(default)s',
                            default=3,
                            type=int)
        parser.add_argument('--tabix', help='Path to tabix executable v0.2.5+.  default: %(default)s',
                            default="tabix")
        parser.add_argument('--chunk_length', help="Chunk length on disk in bytes. default: %(default)s",
                            default=2**19,
                            type=int)
        parser.add_argument('--chunk_width', help="Chunk width on disk in bytes. default: %(default)s",
                            default=64,
                            type=int)
        parser.add_argument('--log', help='Path to logfile, stdout or stderr.  default: stdout',
                            default="stdout")
        parser.add_argument('--zip', action='store_true',
                            help='If flag exists, entire zarr folder is zipped with no compression, ' +
                            "and the original zarr folder is deleted. " +
                            "The zip file name is output_zarr_path appended with '.zip'")
        parser.add_argument('--permissions',
                            help='Unix permissions in octal notation to recursively set on Zarr directory.  ' +
                            'Directories within the Zarr will always be executable to allow traversal.  default: %(default)s',
                            default=0o660,
                            type=auto_int)

        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit()

        args = parser.parse_args()
        input_vcf = args.input_vcf
        output_zarr_path = args.output_zarr_path

        if not input_vcf or not output_zarr_path:
            print ("Missing required options --input_vcf or --output_zarr_path.\n")
            parser.print_help()
            sys.exit()

        group = args.group
        compress_algo = args.compress_algo
        compress_level = args.compress_level
        compress_shuffle = args.compress_shuffle
        alt_number = args.alt_number
        tabix = args.tabix
        chunk_length = args.chunk_length
        chunk_width = args.chunk_width
        is_zip = args.zip
        permissions = args.permissions

        if args.fields:
            fields = [field.strip() for field in args.fields.split(",")]
        else:
            fields = []

        if args.region and "," in args.region:
            region = [single_region.strip() for single_region in args.region.split(",")]
        else:
            region = args.region

        log = args.log
        if args.log:
            if args.log.strip() == "stderr":
                fh_log = sys.stderr
            elif args.log.strip() == "stdout":
                fh_log = sys.stdout
            else:
                fh_log =  open(log, 'w')
        else:
            fh_log = sys.stdout


    with fh_log as handle:
        build_sample_zarr(input_vcf_path=input_vcf,
                        output_zarr_path=output_zarr_path,
                        group=group,
                        region=region,
                        fields=fields,
                        alt_number=alt_number,
                        compress_algo=compress_algo,
                        compress_level=compress_level,
                        compress_shuffle=compress_shuffle,
                        tabix=tabix,
                        chunk_length=chunk_length,
                        chunk_width=chunk_width,
                        log=handle,
                        permissions=permissions)

        if is_zip:
            zip_zarr(zarr_path=output_zarr_path, del_orig=True)
