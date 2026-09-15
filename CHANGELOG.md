# v0.4.3

* Resolved #52 where GuideMaker would crash while calulating Doench efficiency scores if there was an ambiguous nucleotide anywhere in the guide sequence.
    The fix from v0.4.2 was extended to reject guides containing these additional nucleotides: `[WSKMYRVHDB]`.
* Fixed biopython deprecation warnings by updating `feature.strand` to `feature.location.strand`
* Fixed `.gitignore` as it was previously ignoring guidemaker/* due to mangling at the bottom.

# v0.4.2

Fixed a bug where calulating Doench efficiency scores raised an error if there was an 'N' in the first three nucleotides past the PAM in the flanking genomic sequence.  Guidemaker now removes those guides from consideration and reports it as a warning if the flag `--doench_efficiency_score` is used.

# v0.4.1

* Changed how Guidemaker handles DNA sequences that are soft-masked with lowercase letters. The new behavior unmasks all
    sequences and finds guide candidates (and filters them for distance) against the entire sequence.  It also fixes a key value error in computing the Doench efficiency scores when mixed case sequences were used. Previously, guides were only identified in regions where the PAM was all capitalized and edit distance was different if the case was different (a capitalized and lowercase guide were not considered the same).
* Dockerfiles now have more version information
* GitHub build actions were improved

# v0.4.0

*  changed the cli flag `--filter_by_locus` to `--filter_by_attribute`  and added the flag `--attribute_key` so that keys other than the addribute key "locus_tag" can be filtered.
*  Changed the cli to add `--raw_output_only`. This option allows the user to just get the guides that meet LSR and distance criteria without doing any parsing of the genome annotation files.
*  updated caching for Streamlit 1.26
*  updated GC and position methods for Biopython 1.8.1
*  replaced append methods with concat methods for Pandas 2.1.1
*  output data is now gzipped
*  updated Dockerfile to use Minimamba base image
*  Updates to Python dependencies

