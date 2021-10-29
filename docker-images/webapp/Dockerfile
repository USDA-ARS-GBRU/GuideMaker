FROM continuumio/miniconda3:4.9.2
LABEL maintainer="Adam R. Rivers"
LABEL version="0.2.0"
LABEL Description="GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems"

RUN apt-get --allow-releaseinfo-change -y update && apt-get install -y \
    build-essential \
    zlib1g-dev

SHELL ["/bin/bash", "--login", "-c"]

# installation specific to GuideMaker
RUN conda install --yes --channel conda-forge --channel bioconda --channel gagan3012 --freeze-installed \
    nomkl \
    python=3.7 \
    pandas=1.0.3 \
    pybedtools \
    biopython=1.79 \
    numpy=1.19.4 \
    pyyaml=5.4.1 \
    && conda clean -afy

RUN pip install --no-binary :all: nmslib

ENV PATH="~/miniconda3/bin:$PATH"
RUN git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git && cd GuideMaker && pip install -e .

# port on the docker container
EXPOSE 8501 

ENTRYPOINT ["conda", "run", "-n", "base"]

RUN cp /GuideMaker/guidemaker/data/Pseudomonas_aeruginosa.gbk.gz .
RUN cp /GuideMaker/guidemaker/data/Carsonella_ruddii.gbk.gz .

# run web app
CMD streamlit run /GuideMaker/guidemaker/data/app.py --server.maxUploadSize 500