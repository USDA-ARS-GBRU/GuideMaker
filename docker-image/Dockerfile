FROM continuumio/miniconda3:4.8.2
LABEL maintainer="Adam R. Rivers"
LABEL version="0.2.0"
LABEL Description="GuideFinder: globally design guide RNAs for any CRISPR-Cas system in any small genome"

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev
RUN conda create -n guidefinderenv python=3.7 pandas=1.0.3

SHELL ["conda", "run", "-n", "guidefinderenv", "/bin/bash", "-c"]

#RUN conda update -n base -c defaults conda
RUN conda install --channel bioconda pybedtools
RUN conda install -c conda-forge biopython=1.76
RUN pip install --no-binary :all: nmslib

ENV PATH="~/miniconda3/bin:$PATH"

RUN git clone https://github.com/USDA-ARS-GBRU/GuideFinder.git && cd GuideFinder && pip install -e .
ENTRYPOINT ["conda", "run", "-n", "guidefinderenv", "guidefinder"]
CMD ["-h"]
