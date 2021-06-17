---
title: Command Line
category: Access
order: 3
---


### Installation

GuideMaker can be installed from:

3.1. Bioconda: (preferred method because it handles dependencies):

```bash

conda install -c bioconda guidemaker

```

3.2. Github: https://github.com/USDA-ARS-GBRU/GuideMaker

```{bash}
    # Create a conda environment and install and pybedtools
    conda create -n gmenv python=3.7 pybedtools=0.8.2
    conda activate gmenv

    git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git
    cd GuideMaker
    pip install .

    # check if the installation works
    guidemaker -h
```

3.3. Docker image: Available at [Github Registry](https://github.com/orgs/USDA-ARS-GBRU/packages?repo_name=GuideMaker)

```bash

docker pull ghcr.io/usda-ars-gbru/guidemaker-nonavx:sha-9be9fe1c9dca

```


## Dependencies

*   ``pybedtools``
*   ``NMSLib``
*   ``Biopython``
*   ``Pandas``
*   ``Streamlit for webapp``
*   ``altair for plotting``



## Command Line Usage

```{bash}
GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems.

optional arguments:
  -h, --help            show this help message and exit
  --genbank GENBANK [GENBANK ...], -i GENBANK [GENBANK ...]
                        One or more genbank .gbk or gzipped .gbk files for a
                        single genome
  --pamseq PAMSEQ, -p PAMSEQ
                        A short PAM motif to search for, it may use IUPAC
                        ambiguous alphabet
  --outdir OUTDIR, -o OUTDIR
                        The directory for data output
  --pam_orientation {5prime,3prime}, -r {5prime,3prime}
                        PAM position relative to target: 5prime:
                        [PAM][target], 3prime: [target][PAM]. For example,
                        Cas9 is 3prime. Default: '5prime'.
  --guidelength [10-27], -l [10-27]
                        Length of the guide sequence. Default: 20.
  --lsr [0-27]          Length of a seed region near the PAM site required to
                        be unique. Default: 10.
  --dist [0-5]          Minimum hamming distance from any other potential
                        guide. Default: 2.
  --before [1-500]      keep guides this far in front of a feature. Default:
                        100.
  --into [1-500]        keep guides this far inside (past the start site)of a
                        feature. Default: 200.
  --knum [2-20]         how many sequences similar to the guide to report.
                        Default: 3.
  --controls CONTROLS   Number or random control RNAs to generate. Default:
                        1000.
  --threads THREADS     The number of cpu threads to use. Default: 2
  --log LOG             Log file
  --tempdir TEMPDIR     The temp file directory
  --restriction_enzyme_list [RESTRICTION_ENZYME_LIST [RESTRICTION_ENZYME_LIST ...]]
                        List of sequence representing restriction enzymes.
                        Default: None.
  --keeptemp            Option to keep intermediate files be kept
  --plot                Option to genereate guidemaker plots
  --config CONFIG       Path to YAML formatted configuration file, default is 
                        /Users/admin/opt/anaconda3/envs/gmenv/lib/python3.7/si
                        te-packages/guidemaker/data/config_default.yaml
  -V, --version         show program's version number and exit

To run the web app locally, in terminal run:
-----------------------------------------------------------------------
streamlit run /Users/admin/opt/anaconda3/envs/gmenv/lib/python3.7/site-
packages/guidemaker/data/app.py
-----------------------------------------------------------------------

```

## Examples

Use case: Make 20 nucleotide guide sequences for SpCas9 (NGG) in the bacterium
__Carsonela ruddii__. The length of the seed region near the PAM required to be
unique in each guide is 11 nucleotides.

```{bash}
    guidemaker \
    -i tests/test_data/Carsonella_ruddii.gbk \
    -p NGG \
    --pam_orientation 3prime \
    --guidelength 20 \
    --lsr 11 \
    -o OUTDIR \
    --threads 2

```
> [Questions](http://example.com/signup) or suggestions- Let us know at [example.com](http://example.com/).

![](//placehold.it/800x600)
