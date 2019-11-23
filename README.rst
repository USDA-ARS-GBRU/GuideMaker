predictPAM: A python module to predict custom PAM sites in any small genome
==================================================================================================

Authors
-------
* Ravin Poudel, US Department of Agriculture, Agricultural Research Service
* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service
* Christopher Reisch, Department of Microbiology and Cell Science, University of Florida


Introduction
-------------



Installation
-------------
predictPAM can be installed from:

1. The Github repository: https://github.com/USDA-ARS-GBRU/predictPAM

.. code-block:: bash

    git clone https://github.com/USDA-ARS-GBRU/predictPAM.git


Dependencies
-------------
Following are the required softwares/programs.

- ``seqkit``

- ``pybedtools``

- ``NMSLIB``

- ``Biopython``


Usage
---------

-h, --help            	Show this help message and exit.

--gbkfile, -i		    A ``.GBk``, ``.GB`` file or files. Supports single or multiple genome files with single or multiple chromosomes. Required.

--pamseq, -p			A short PAM motif to search for, may be use IUPAC ambiguous alphabet. Required.

--targetlength, -l      Length of the target sequence. Default is 22.

--strand, -s            Strand of DNA to search for PAM motif. Default is forward. Options: {forward, reverse}.

--lcp                   Length of conserved sequence close to PAM motif. Default is 12.

--eds                   Unexcepted Levenshtein edit distance on the distal portion of target sequence from PAM motif. Options: {0, 1, 2, 3, 4, 5}. Default is 2.

--outfile, -o           The table of pam sites and data. Required.

--tempdir				Specify the temp file directory. Default is None.

--keeptemp				Should intermediate files be kept? Default is false.

--log		          	Log file. Default is predictPAM.log.

--threads		     	Number of processor threads to use. Default is 1.


Examples
---------

Use case 1: Retrieving target sequence for a given PAM motif in the forward strand, where length of target sequence is 25 base pair.
12 base pair close to PAM motif is conserved ,i.e. unique. Remaining sequence at the distal end to PAM motif has leven distance of more than 2.
Return a table of pam sites and associated data.

.. code-block:: bash
    
    predictPAM -i sample.gbk -p ATCGAT --targetlength 25 -strand forward \
    --lcp 12 --eds 2 --outfile out.txt \
    --log logfile.txt --keeptemp -threads 2

ITSxpress can take gzipped or un-gzipped FASTQ files and it can write gzipped or
un-gzipped FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz or fastq.gz.

Use case 2: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
forward and reverse gzipped FASTQ files using two cpu threads. Return a forward
and reverse read files  for use in Dada2.

.. code-block:: bash

    itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz --region ITS2 \
    --taxa Fungi --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

ITSxpress can take gzipped or un-gzipped FASTQ files and it can write gzipped or
un-gzipped FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz or fastq.gz.


Use case 3: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
an interleaved gzipped FASTQ files using two cpu threads. Return a single merged file for use in Deblur.

.. code-block:: bash

    itsxpress --fastq interleaved.fastq.gz  --region ITS2 --taxa Fungi \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2


Use case 4: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
an single-ended gzipped FASTQ files using two cpu threads.

.. code-block:: bash

    itsxpress --fastq single-end.fastq.gz --single_end --region ITS2 --taxa Fungi \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

Single ended data is less common and may come from a dataset where the reads have already
been merged.

Use case 5: Trimming the ITS1 region from a Alveolata amplicon sequencing dataset with
an interleaved gzipped FASTQ files using 8 cpu threads.

.. code-block:: bash

    itsxpress --fastq interleaved.fastq.gz --region ITS1 --taxa Alveolata \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 8


License information
--------------------
This software is a work of the United States Department of Agriculture,
Agricultural Research Service and is released under a Creative Commons CC0
public domain attribution.
