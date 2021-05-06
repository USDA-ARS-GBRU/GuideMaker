FROM continuumio/miniconda3:4.9.2
LABEL maintainer="Adam R. Rivers"
LABEL version="0.2.0"
LABEL Description="GuideMaker: Globally design guide RNAs for any CRISPR-Cas system in any small genome"

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev

SHELL ["/bin/bash", "--login", "-c"]

# installation specific to GuideMaker
RUN conda install --yes --channel conda-forge --channel bioconda --channel gagan3012 --freeze-installed \
    nomkl \
    python=3.7 \
    pybedtools=0.8.2 \
    && conda clean -afy

RUN pip install --no-binary :all: nmslib

ENV PATH="~/miniconda3/bin:$PATH"
RUN git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git && cd GuideMaker && pip install -e .

EXPOSE 8501

ENTRYPOINT ["conda", "run", "-n", "base"]

#CMD streamlit run ~/GuideMaker/guidemaker/data/app.py

#CMD ls /GuideMaker/guidemaker/data/*


CMD streamlit run /GuideMaker/guidemaker/data/app.py


#CMD streamlit run --server.port 8501 --server.enableCORS false /GuideMaker/guidemaker/data/app.py


#CMD conda list 

# docker run -p 8501:8501 gmapp