#GuideMaker
## Globally design gRNAs for any CRISPR-Cas system in any small genome


## Authors

* Lidimarie Trujillo, Department of Microbiology and Cell Science, University of Florida
* Ravin Poudel, PhD, Department of Microbiology and Cell Science, University of Florida
* Christopher Reisch, PhD, Department of Microbiology and Cell Science, University of Florida
* Adam R. Rivers, PhD , US Department of Agriculture, Agricultural Research Service


## Introduction



## Installation

GuideMaker can be installed from:

1. The Github repository: https://github.com/USDA-ARS-GBRU/guidemaker

```{bash}
    git clone https://github.com/USDA-ARS-GBRU/GuideMaker.git
```


## Dependencies

Following are the required softwares/programs.

* ``pybedtools``
* ``NMSLib``
* ``Biopython``
* ``Pandas``


## Usage

```
GuideMaker: globally design guide RNAs for any CRISPR-Cas system in any genome

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
  --strand {forward,reverse,both}, -s {forward,reverse,both}
                        Strand of DNA to search
  --lcp [0-27]          Length of the guide closest to the PAM required to be
                        unique
  --dist [1-5]          Minimum hamming distance from any other potential
                        guide
  --before [1-500]      keep guides this far in front of a feature
  --into [1-500]        keep guides this far inside (past the start site)of a
                        feature
  --knum [2-20]         how many sequences similar to the guide to report
  --controls CONTROLS   The number or random control RNAs to generate
  --threads THREADS     The number of cpu threads to use
  --log LOG             Log file
  --tempdir TEMPDIR     The temp file directory
  -V, --version         show program's version number and exit
```

## Examples


Use case: Retrieving target sequence for a given PAM motif in the forward and reverse strands, where length of guide sequence is 20 base pair.
12 base pair close to PAM motif is conserved ,i.e. unique and the full sequence has a hamming distance of more than 2.
Here the number of used threads is 2
Return a table of pam sites and associated data, at in current folder.

```
		guidemaker -i test/test_data/Burkholderia_thailandensis_E264_ATCC_700388_133.gbk \
		--pamseq AGG  --outdir OUTDIR --pam_orientation 5prime \
		--guidelength 20 --strand both --lcp 10 --dist 3 --before 100 \
		--into  100 --knum 10 --controls 10  --log logfile.txt --threads 2

```

## API documentation

API documentation for the module can be found [here](/html/guidemaker/index.html)

## License information

As a work of the United State Governemnt this software is available under  CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
