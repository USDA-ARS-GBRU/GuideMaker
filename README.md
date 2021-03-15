
# GuideMaker: Globally design guide RNAs for any CRISPR-Cas system in any small genome


## Authors


* Ravin Poudel, PhD, Department of Microbiology and Cell Science, University of Florida
* Lidimarie Trujillo, Department of Microbiology and Cell Science, University of Florida
* Christopher Reisch, PhD, Department of Microbiology and Cell Science, University of Florida
* Adam R. Rivers, PhD , US Department of Agriculture, Agricultural Research Service



## Installation

GuideMaker can be installed from:

1. The Github repository: https://github.com/USDA-ARS-GBRU/GuideMaker

```{bash}
    git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git
```


## Dependencies

* ``pybedtools``
* ``NMSLib``
* ``Biopython``
* ``Pandas``
* ``Streamlit for webapp``
* ``altair for plotting``


## Usage

```
GuideMaker: Globally design guide RNAs for any CRISPR-Cas system in any small genome

optional arguments:

  -h, --help            show this help message and exit

  --genbank GENBANK [GENBANK ...], -i GENBANK [GENBANK ...]
                        One or more genbank .gbk or gzipped .gbk files for a
                        single genome

  --pamseq PAMSEQ, -p PAMSEQ
                        A short PAM motif to search for, it may use IUPAC
                        ambiguous alphabet
  --outfile OUTFILE, -o OUTFILE
                        The table of PAM sites and data
  --pam_orientation {5prime,3prime}, -r {5prime,3prime}
                        PAM position relative to target: 5prime:
                        [PAM][target], 3prime: [target][PAM]. For example,
                        Cas9 is 3prime
  --guidelength [10-27], -l [10-27]
                        Length of the guide sequence
  --lsr [0-27]          Length of a seed region near the PAM site required to
                        be unique
  --dist [0-5]          Minimum hamming distance from any other potential
                        guide. Default dist is >= 2.
  --before [1-500]      keep guides this far in front of a feature
  --into [1-500]        keep guides this far inside (past the start site)of a
                        feature
  --knum [2-20]         Number of sequences similar to the guide to report
  --controls CONTROLS   Number or random control RNAs to generate
  --threads THREADS     Number of cpu threads to use
  --log LOG             Log file
  --tempdir TEMPDIR     The temp file directory
  --restriction_enzyme_list List of sequence representing restriction enzymes
  --keeptemp            Option to keep intermediate files
  --plot                Option to generate plots
  --config CONFIG       Path to YAML formatted configuration file
  -V, --version         show program's version number and exit
```

## Examples


Use case: Make 20 nucleotide guide sequences for SpCas9 (NGG) in Carsonela ruddii. The length of the seed region near the PAM requred to be unique in each guide is 11 nucleotides.

```
    guidemaker \
    -i test/test_data/Carsonella_ruddii.gbk \
    -p NGG \
    --pam_orientation 3prime \
    --guidelength 20 \
    --lsr 11 \
    -o OUTDIR \
    --threads 2

```

### Running Web App locally

```{bash}
streamlit run guidemaker/app.py 
```

[![IMAGE ALT TEXT HERE](https://github.com/USDA-ARS-GBRU/GuideMaker/blob/rp01/GuideMakerApp.png)](https://guidemaker.org)

## API documentation

API documentation for the module can be found [here](https://guidemaker.org/html/guidemaker/index.html)

## License information

Guidemaker was created by the [United States Department of Agriculture - Agricultural Research Service (USDA-ARS)](https://www.ars.usda.gov/) this software is available under [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/)

