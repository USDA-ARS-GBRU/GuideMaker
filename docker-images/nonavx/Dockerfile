FROM continuumio/miniconda3:4.9.2
LABEL maintainer="Adam R. Rivers"
LABEL version="0.2.0"
LABEL Description="GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems"

RUN apt-get --allow-releaseinfo-change -y update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

SHELL ["/bin/bash", "--login", "-c"]

RUN conda install --yes --channel conda-forge --channel bioconda --freeze-installed \
    nomkl \
    python=3.7 \
    pandas=1.0.3 \
    pybedtools \
    biopython=1.79 \
    numpy=1.19.4 \
    pyyaml=5.4.1 \
    nmslib \
    jsonschema==3.2 \
    && conda clean -afy

ENV PATH="~/miniconda3/bin:$PATH"
RUN git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git && cd GuideMaker && pip install .
RUN mkdir /data
WORKDIR /data
ENTRYPOINT ["conda", "run", "-n", "base", "guidemaker"]


