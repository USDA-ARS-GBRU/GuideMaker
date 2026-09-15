FROM ghcr.io/mamba-org/micromamba:1.5.1-bullseye AS guidemaker
LABEL maintainer="Adam R. Rivers"
LABEL Description="GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems"

RUN micromamba install --channel conda-forge --channel bioconda -y -n base \
    nomkl \
    python=3.9 \
    pip=23.2\
    pybedtools=0.9.1 \
    && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential=12.9 \
    python3-dev=3.9.2-3 \
    zlib1g-dev=1:1.2.11.dfsg-2+deb11u2 \
    git=1:2.30.2-1+deb11u2 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

RUN pip install --no-binary :all: nmslib

RUN git clone https://github.com/mjfos2r/GuideMaker.git && cd GuideMaker && pip install .

# port on the docker container
EXPOSE 8501 


ENTRYPOINT ["micromamba", "run", "-n", "base"]

RUN cp GuideMaker/guidemaker/data/Pseudomonas_aeruginosa.gbk.gz .
RUN cp GuideMaker/guidemaker/data/Carsonella_ruddii.gbk.gz .

# run web app
CMD streamlit run GuideMaker/guidemaker/data/app.py --server.maxUploadSize 5000
