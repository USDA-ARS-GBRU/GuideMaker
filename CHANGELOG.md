# v0.4.0

*  changed the cli flag `--filter_by_locus` to `--filter_by_attribute`  and added the flag `--attribute_key` so that keys other than the addribute key "locus_tag" can be filtered.
*  Changed the cli to add `--raw_output_only`. This option allows the user to just get hte guides that met LSR and distance criteria without doing any parsing of the genome annotation files.
*  updated caching for Streamlit 1.26
*  updated GC and position methods for Biopython 1.8.1
*  replaced append methods with concat methods for Pandas 2.1.1
*  output data is now gzipped
*  updated Dockerfile to use minimamba base image
*  Updates to python dependencies

