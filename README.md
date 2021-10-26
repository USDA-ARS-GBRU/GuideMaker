<h1 style="color: #FA4616" >GuideMaker</h1>

**GuideMaker: Software to design CRISPR-Cas guide RNA pools in non-model genomes** 🦠 🧬

CRISPR-Cas systems have expanded the possibilities for gene editing in bacteria and eukaryotes. There are many excellent tools for designing the CRISPR-Cas guide RNAs for model organisms with standard Cas enzymes. GuideMaker is intended as a fast and easy-to-use design tool for atypical projects with 1) non-standard Cas enzymes, 2) non-model organisms, or 3) projects that need to design a panel of guide RNAs (gRNA) for genome-wide screens.

GuideMaker can rapidly design gRNAs for gene targets across the genome from a degenerate protospacer adjacent motif (PAM) and a GenBank file or Fasta and GFF/GTF file. The tool applies Hierarchical Navigable Small World (HNSW) graphs to speed up the comparison of guide RNAs enabling the user to design gRNAs for all genes for a typical bacterial genome and PAM sequence in about 1-2 minutes on a laptop.

GuideMaker enables the rapid design of genome-wide CRISPR/Cas gene function studies in non-model organisms with any Cas enzyme. While GuideMaker is designed with prokaryotic genomes in mind, it can process smaller eukaryotic genomes as well. GuideMaker is available as command-line software and as a **[web application](https://guidemaker.app.scinet.usda.gov)** at **https://guidemaker.app.scinet.usda.gov** and in the **[CyCverse Discovery Environment](https://cyverse.org/discovery-environment)**.

## Methods to access GuideMaker

GuideMaker can be easily accessed via:
- Web Application
- CyVerse Discovery Environment
- Command Line 
- Local Web Application

**NOTE:** *Our web application runs on a small server instance and is primarily designed for the lower-memory requirements bacterial genomes. We recommend that users run larger genomes on the **[CyCverse Discovery Environment](https://de.cyverse.org/apps/de/518589c0-994a-11ea-9ea3-008cfa5ae621)**  or run GuideMaker locally as a command-line or web browser-based application.

1.[Web Application](https://guidemaker.app.scinet.usda.gov)|  2.[CyCverse Discovery Environment](https://de.cyverse.org/apps/de/518589c0-994a-11ea-9ea3-008cfa5ae621)
:-------------------------:|:-------------------------:
[![Image of GuideMaker Web App](https://raw.githubusercontent.com/USDA-ARS-GBRU/GuideMaker/main/guidemaker/data/scinet.png)](https://guidemaker.app.scinet.usda.gov)|[![Image of GuideMaker Web App](https://raw.githubusercontent.com/USDA-ARS-GBRU/GuideMaker/main/guidemaker/data/cyverse.png)](https://de.cyverse.org/apps/de/518589c0-994a-11ea-9ea3-008cfa5ae621)


## 3.Command Line

GuideMaker can be installed from:

### 3.1. Bioconda: (preferred method because it handles dependencies):

```bash
# Create a conda environment and install GuideMaker via Bioconda.

conda create --strict-channel-priority --override-channels --channel conda-forge --channel bioconda --channel defaults --name gmenv guidemaker

# Activate conda env
conda activate gmenv

# Test the installation
guidemaker -h

```

### 3.2. [Github](https://github.com/USDA-ARS-GBRU/GuideMaker)

```bash
    # Create a conda environment and install and pybedtools
    conda create -n gmenv python=3.7 pybedtools=0.8.2
    conda activate gmenv

    git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git
    cd GuideMaker
    pip install .

    # check if the installation works
    guidemaker -h
```

### 3.3. Docker image: Available at [Github Registry](https://github.com/orgs/USDA-ARS-GBRU/packages?repo_name=GuideMaker)

```bash

docker pull ghcr.io/usda-ars-gbru/guidemaker-nonavx:sha-9be9fe1c9dca

```


### Dependencies

*   ``pybedtools``
*   ``NMSLib``
*   ``Biopython``
*   ``Pandas``
*   ``Streamlit for webapp``
*   ``altair for plotting``



### Command Line Usage

```
usage: guidemaker [-h] --genbank GENBANK [GENBANK ...] --pamseq PAMSEQ
                  --outdir OUTDIR [--pam_orientation {5prime,3prime}]
                  [--guidelength [10-27]] [--lsr [0-27]]
                  [--dtype {hamming,leven}] [--dist [0-5]] [--before [1-500]]
                  [--into [1-500]] [--knum [2-20]] [--controls CONTROLS]
                  [--threads THREADS] [--log LOG] [--tempdir TEMPDIR]
                  [--restriction_enzyme_list [RESTRICTION_ENZYME_LIST [RESTRICTION_ENZYME_LIST ...]]]
                  [--filter_by_locus [FILTER_BY_LOCUS [FILTER_BY_LOCUS ...]]]
                  [--doench_efficiency_score] [--cfd_score] [--keeptemp]
                  [--plot] [--config CONFIG] [-V]

GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas
systems

optional arguments:
  -h, --help            show this help message and exit
  --genbank GENBANK [GENBANK ...], -i GENBANK [GENBANK ...]
                        One or more genbank .gbk or gzipped .gbk files for a single genome. Provide this or GFF and fasta files
  --fasta FASTA [FASTA ...], -f FASTA [FASTA ...]
                        One or more fasta or gzipped fasta files for a single genome. If using a fasta, a GFF must also be provided but not a genbank file.
  --gff GFF [GFF ...], -g GFF [GFF ...]
                        One or more genbank GFF files for a single genome. If using a GFF a fasta file must also be provided but not a genbank file.
  --pamseq PAMSEQ, -p PAMSEQ
                        A short PAM motif to search for, it may use IUPAC ambiguous alphabet
  --outdir OUTDIR, -o OUTDIR
                        The directory for data output
  --pam_orientation {5prime,3prime}, -r {5prime,3prime}
                        PAM position relative to target: 5prime: [PAM][target], 3prime: [target][PAM]. For example, Cas9 is 3prime. Default: '5prime'.
  --guidelength [10-27], -l [10-27]
                        Length of the guide sequence. Default: 20.
  --lsr [0-27]          Length of a seed region near the PAM site required to be unique. Default: 10.
  --dtype {hamming,leven}
                        Select the distance type. Default: hamming.
  --dist [0-5]          Minimum edit distance from any other potential guide. Default: 2.
  --before [1-500]      keep guides this far in front of a feature. Default: 100.
  --into [1-500]        keep guides this far inside (past the start site)of a feature. Default: 200.
  --knum [2-20]         how many sequences similar to the guide to report. Default: 5.
  --controls CONTROLS   Number or random control RNAs to generate. Default: 1000.
  --threads THREADS     The number of cpu threads to use. Default: 2
  --log LOG             Log file
  --tempdir TEMPDIR     The temp file directory
  --restriction_enzyme_list [RESTRICTION_ENZYME_LIST [RESTRICTION_ENZYME_LIST ...]]
                        List of sequence representing restriction enzymes. Default: None.
  --filter_by_locus [FILTER_BY_LOCUS [FILTER_BY_LOCUS ...]]
                        List of locus tag. Default: None.
  --doench_efficiency_score
                        Doench et al. 2016 - only for NGG PAM: Default: None.
  --cfd_score           CFD score for assessing off-target activity of gRNAs: Default: None.
  --keeptemp            Option to keep intermediate files be kept
  --plot                Option to genereate guidemaker plots
  --config CONFIG       Path to YAML formatted configuration file, default is /Users/rivers/Documents/guidemaker/guidemaker/data/config_default.yaml
  -V, --version         show program's version number and exit


To run the web app locally, in terminal run:
-----------------------------------------------------------------------
streamlit run [Dynamically created path to app file]
-----------------------------------------------------------------------


```

### Examples

Use case: Make 20 nucleotide guide sequences for SpCas9 (NGG) in the bacterium
_Carsonela ruddii_. The length of the seed region near the PAM required to be
unique in each guide is 11 nucleotides.

```bash
    guidemaker \
    -i tests/test_data/Carsonella_ruddii.gbk \
    -p NGG \
    --pam_orientation 3prime \
    --guidelength 20 \
    --lsr 11 \
    -o OUTDIR \
    --doench_efficiency_score \
    --threads 2 

```

## 4. Running Web App locally
To run the web app locally, you first need to complete the command line installation described above.

If the path of the `app.py` differs from the one displayed below, you can locate the path by first running `guidemaker --help`. Script to run the web app locally is available at the bottom of the help command output. 

```bash

streamlit run /[user path prefix]/anaconda3/envs/gmenv/lib/python3.7/site-packages/guidemaker/data/app.py --server.maxUploadSize 500

```
![Image of GuideMaker Web App](https://raw.githubusercontent.com/USDA-ARS-GBRU/GuideMaker/main/guidemaker/data/scinet.png)

## Using GuideMaker's results

_This section provides information on how to use GuideMaker's results to create a molecular protocol for pooled CRISPR screens._

<h4
	style="color: #003087"
	>Pooled CRISPR Experiments
</h4>

Experiments that target the entire genome, or many genes at once, are typically performed in pooled experiments where 100-100,000+ targets are tested simultaneously. The pooled oligonucleotides for each gRNA are cloned in one batch and used simultaneously in the designed experiment. Each gRNA sequence acts as a barcode that can be quantified with high-throughput sequencing to elucidate each target's relative importance under the experimental conditions. 


<h4 style="color: #003087" >Vectors for gRNA Cloning</h4>

Genome-scale CRISPR experiments require a gRNA vector amenable to high-throughput cloning, most often through [Golden Gate cloning](https://blog.addgene.org/plasmids-101-golden-gate-cloning), a restriction enzyme-dependent reaction. Plasmids to express gRNA are available from Addgene can be found at the link below, though not all of these are compatible with high-throughput cloning. 


<h4 style="color: #003087" >Addgene: CRISPR Plasmids - Empty gRNA Vectors</h4>

After running GuideMaker, the designed gRNA output can be downloaded and with minor adjustments, the targets can be ordered as oligos for cloning. Pooled oligonucleotides can be purchased from several vendors, including those listed below. Pool sizes vary from 100 to over 200,000 oligonucleotides. Vendor specifications for the number of oligos, oligo length, and cost per bp vary widely. For bacterial genome-scale experiments, as of 2021, Genscript offers pool sizes of 12,472 and 91,766 with up to 79 bp per oligo for list prices of $1600 and $4,000, respectively.

Some example vendors are:

*   [GenScript](https://www.genscript.com/precise-synthetic-oligo-pools.html)
*   [Twist Bioscience](https://www.twistbioscience.com/products/oligopools)
*   [Agilent](https://www.agilent.com/en/product/sureprint-oligonucleotide-library-synthesis/oligonucleotide-library-synthesis/sureprint-oligonucleotide-libraries-288039)
*   [Arbor Biosciences](https://arborbiosci.com/oligos-and-arrays/dna-and-rna-oligo-pools/)

Most pools require amplification before cloning to convert the ssDNA to dsDNA and increase the concentration for efficient cloning. Accordingly, adding a constant region at the 3' end for primer binding is recommended. Sub-pools can also be amplified by adding unique constant regions to some oligos, enabling the large-scale synthesis to be split amongst organisms or specific targets in a single organism. Because Golden Gate cloning utilizes restriction enzymes, filtering gRNA designs with the cognate restriction enzyme recognition sites is necessary, a feature found in GuideMaker. A general protocol for cloning pooled gRNA from synthesized oligonucleotides from IDT is linked below, though similar workflows can be used for pools from other vendors. 

*   [Cloning high-quality CRISPR libraries with oPool Oligo Pools (SYB-10182-PR12/2019)](https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/user-submitted-method/cloning-high-quality-crispr-libraries-with-opools-oligo-pools-user-method.pdf?sfvrsn=3db31607_7)
*   [Addgene: Guide to Using Pooled Libraries](https://www.addgene.org/guides/pooled-libraries/)

<h4 style="color: #003087">Pooled CRISPR Data Analysis</h4>
After the experiment, the cells are collected and DNA is isolated. The target sequence is then amplified and adaptors for high-throughput sequencing added. Several data analysis pipelines have been developed to identify target sequences over-represented or under-represented in the pool. The manuscript by Wang et al. (2019) provides a protocol for using a high-quality tool with these capabilities. 

<h4 style="color: #003087">Citation</h4>
Wang, B., Wang, M., Zhang, W. et al. Integrative analysis of pooled CRISPR genetic screens using MAGeCKFlute. Nat Protoc 14, 756–780 (2019). https://doi.org/10.1038/s41596-018-0113-7

### FAQs

Coming soon...

### Reporting Errors and Suggestions
Open the GuideMaker [github repo](https://github.com/USDA-ARS-GBRU/GuideMaker), navigate to the `Issues` page and submit an `issue` to report difficulties, errors, or suggestions for improvements. Also, check **FAQs** section prior submitting an issue. 

![Image of GuideMaker Web App](https://raw.githubusercontent.com/USDA-ARS-GBRU/GuideMaker/main/guidemaker/data/gitissue.png)


## Citation
*Poudel R, Rodriguez LT, Reisch CR, Rivers AR. GuideMaker: Software to design CRISPR-Cas guide RNA pools in non-model genomes. 2021. In Review.

## API documentation

API documentation for the module can be found [here](https://ravinpoudel.github.io/GuideMaker/index.html)

## License information
GuideMaker was created by the [United States Department of Agriculture - Agricultural Research Service 
(USDA-ARS)](https://www.ars.usda.gov/). As a work of the United States Government, this software is available under 
the [CC0 1.0 Universal Public Domain Dedication (CC0 1.0)](https://creativecommons.org/publicdomain/zero/1.0)


[![CircleCI](https://img.shields.io/circleci/build/github/USDA-ARS-GBRU/GuideMaker?logo=CircleCi&token=802d114b3ec676d153b4b9fa6a781f9345756fc9)](https://app.circleci.com/pipelines/github/USDA-ARS-GBRU/GuideMaker)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/0f49664d414e44159c1f195474027eae)](https://www.codacy.com/gh/USDA-ARS-GBRU/GuideMaker/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=USDA-ARS-GBRU/GuideMaker&amp;utm_campaign=Badge_Grade)
[![Codecov](https://img.shields.io/codecov/c/github/USDA-ARS-GBRU/GuideMaker?logo=codecov)](https://app.codecov.io/gh/USDA-ARS-GBRU/GuideMaker)
[![DOI](https://zenodo.org/badge/217529920.svg)](https://zenodo.org/badge/latestdoi/217529920)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/guidemaker/badges/downloads.svg)](https://anaconda.org/bioconda/guidemaker)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/USDA-ARS-GBRU/GuideMaker?style=social)

## About us

GuideMaker was developed by the USDA Agricultural Research Service, Genomics and Bioinformatics Research Unit group in Gainesville, FL led by Adam Rivers. Check out our other work at https://tinyecology.com.
