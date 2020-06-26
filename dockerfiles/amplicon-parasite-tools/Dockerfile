FROM continuumio/miniconda3:4.8.2

# Create a yml file specifying a conda environment that contains all tools for the amplicon SNP-calling (BCL-to-VCF, "Stage1, Steps 1-3") parasite pipeline
# (except for those tools used to process manifest files, see https://gitlab.com/malariagen-aspis/aspis-pipeline/-/blob/e8b35283ad70c41b4d0c9c9a3a660e31fff4431b/Docker/ManifestTools/Dockerfile).
# Copied from https://gitlab.com/malariagen-aspis/aspis-pipeline/-/blob/e8b35283ad70c41b4d0c9c9a3a660e31fff4431b/Stage1-Pipeline/conda/pipe-tools.txt.
# Could possibly be simplified or updated further.
# Note that using the Docker images for individual tools (located at https://gitlab.com/malariagen-aspis/aspis-pipeline/-/tree/master/Docker) led to errors from bamcollate2;
# Similarly, using a conda environment with the latest versions of packages from the https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/ channel does not give a working version of bcftools.
RUN echo "# This file may be used to create an environment using:" > amplicon-parasite-tools.yml && \
    echo "# $ conda create --name <env> --file <this file>" >> amplicon-parasite-tools.yml && \
    echo "# platform: linux-64" >> amplicon-parasite-tools.yml && \
    echo "@EXPLICIT" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/teepot-1.2.0-1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/libbzip2-1.0.6-1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/libgcc-ng-8.2.0-hdf63c60_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/libgcc_npg-7.3.0-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/liblzma-5.2.3-1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/libstdcxx-ng-8.2.0-hdf63c60_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/nettle-3.3-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/openssl-1.0.2o-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/zlib-1.2.11-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/bwa-aligner-0.7.17-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/bzip2-1.0.6-1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/gmplib-6.1.2-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/libffi-3.2.1-hd88cf55_4.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/liblief-0.9.0-h7725739_2.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/libpng-1.6.34-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/libxml2-2.9.7-2.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/lz4-c-1.8.1.2-h14c3975_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/lzo-2.10-h49e0be7_2.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/ncurses-6.1-he6710b0_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/patchelf-0.9-he6710b0_3.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/tk-8.6.8-hbc83047_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/xz-5.2.3-1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/yaml-0.1.7-had09818_2.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/gnutls-3.4.17-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/libarchive-3.3.2-h4a23435_4.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/libedit-3.1.20181209-hc058e9b_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/libgd-2.2.5-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/readline-7.0-h7b6447c_5.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/libcurl-7.58.0-1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/sqlite-3.27.2-h7b6447c_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/htslib-1.8-plugins_201712_2.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/python-3.6.5-hc3d631a_2.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/staden_io_lib-1.14.9-4.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/asn1crypto-0.24.0-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/bambi-0.11.1-1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/bcftools-1.8-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/certifi-2019.3.9-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/chardet-3.0.4-py36_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/cryptography-vectors-2.3.1-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/filelock-3.0.10-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/glob2-0.6-py36_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/idna-2.8-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/libmaus2-2.0.420-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/markupsafe-1.1.1-py36h7b6447c_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pkginfo-1.5.0.1-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/psutil-5.6.1-py36h7b6447c_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/py-lief-0.9.0-py36h7725739_2.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pycosat-0.6.3-py36h14c3975_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pycparser-2.19-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pysocks-1.6.8-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/python-libarchive-c-2.8-py36_6.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pytz-2018.9-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pyyaml-5.1-py36h7b6447c_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/ruamel_yaml-0.15.46-py36h14c3975_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/samtools-1.8-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/six-1.12.0-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/soupsieve-1.8-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/staden_io_bin-1.14.9-4.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/tqdm-4.31.1-py36_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/beautifulsoup4-4.7.1-py36_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/linux-64/biobambam2-2.0.79-0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/cffi-1.12.2-py36h2e261b9_1.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/setuptools-40.8.0-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/cryptography-2.3-py36hb7f436b_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/jinja2-2.10-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/wheel-0.33.1-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pip-19.0.3-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/pyopenssl-19.0.0-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/urllib3-1.24.1-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/requests-2.21.0-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/conda-4.6.11-py36_0.tar.bz2" >> amplicon-parasite-tools.yml && \
    echo "https://repo.anaconda.com/pkgs/main/linux-64/conda-build-3.17.8-py36_0.tar.bz2" >> amplicon-parasite-tools.yml

RUN conda create -n amplicon-parasite-tools --file amplicon-parasite-tools.yml && \
    conda run -n amplicon-parasite-tools conda list && \
    conda clean -afy && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete

ENV PATH='/opt/conda/envs/amplicon-parasite-tools/bin:/opt/conda/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin' \
    CONDA_PREFIX='/opt/conda/envs/amplicon-parasite-tools' \
    CONDA_DEFAULT_ENV='amplicon-parasite-tools' \
    CONDA_PROMPT_MODIFIER='(amplicon-parasite-tools) '

RUN echo ". /opt/conda/etc/profile.d/conda.sh" > ~/.bashrc && \
    echo "conda activate amplicon-parasite-tools" >> ~/.bashrc
