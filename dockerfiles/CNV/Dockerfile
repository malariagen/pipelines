FROM --platform=linux/amd64 continuumio/miniconda3

# Update the package lists and install necessary packages
RUN apt-get update && \
    apt-get install -y curl bzip2 ca-certificates gnupg && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#TODO: Determine if R is really required here
# # Add the R repository and install the latest R
# RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu hirsute-cran40/" >> /etc/apt/sources.list && \
#     apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
#     apt-get update && \
#     apt-get install -y r-base && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/*

# Copy the environment file to the container
COPY requirements_conda.yml requirements_conda.yml

# Create a conda environment and activate it
RUN conda env create -f requirements_conda.yml && \
    echo "conda activate $(head -1 /opt/conda/envs/*/etc/environment.yml | cut -d'/' -f6)" >> ~/.bashrc

# Set the working directory
WORKDIR /cnv

# Copy the source code to the container
COPY . .

# Set the default command to run when the container starts
# CMD ["python", ""]