# v0.4.1

* Changed how Guidemaker handles DNA sequenses that are soft-masked with lower case letters. The new behavior unmasks all 
    sequences and finds guide candidates (and filters them for distance) against the entire sequence.  It also fixes a a key value error in computing the doench efficency scores when mixed case sequences were used. Prevously, guides were only identified in regions where the PAM was all capital and edit distance was different if the case was different (a capitalized and lowercase guide were not considered the same).
* Dockerfiles now have more version information
* Github build actions were impoved

# v0.4.0

*  changed the cli flag `--filter_by_locus` to `--filter_by_attribute`  and added the flag `--attribute_key` so that keys other than the addribute key "locus_tag" can be filtered.
*  Changed the cli to add `--raw_output_only`. This option allows the user to just get hte guides that met LSR and distance criteria without doing any parsing of the genome annotation files.
*  updated caching for Streamlit 1.26
*  updated GC and position methods for Biopython 1.8.1
*  replaced append methods with concat methods for Pandas 2.1.1
*  output data is now gzipped
*  updated Dockerfile to use minimamba base image
*  Updates to python dependencies

