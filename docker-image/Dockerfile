FROM continuumio/miniconda3:4.8.2
LABEL maintainer="Adam R. Rivers"
LABEL version="0.2.0"
LABEL Description="GuideMaker: Globally design guide RNAs for any CRISPR-Cas system in any small genome"

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev

SHELL ["/bin/bash", "--login", "-c"]

RUN conda install --yes --channel conda-forge --channel bioconda --freeze-installed \
    nomkl \
    python=3.7 \
    pandas=1.0.3 \
    pybedtools \
    biopython=1.76 \
    numpy=1.19.4 \
    pyyaml==5.3.1 \
    && conda clean -afy

RUN pip install --no-binary :all: nmslib

ENV PATH="~/miniconda3/bin:$PATH"
RUN git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git && cd GuideMaker && pip install -e .
RUN mkdir /data
WORKDIR /data
ENTRYPOINT ["conda", "run", "-n", "base", "guidemaker"]


