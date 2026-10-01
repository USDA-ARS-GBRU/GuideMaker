# Objective

This document describes a refactor of step2.rs and lib.2s to update the protocol in guidemaker-step2  analysis run.  Only thse two fiels should be changed.

## Crurent opperation 
1. opens polars dataframe from inital guidemaker-scan
2. dereplicates reads in the LSR seed region
3.  finds the intervals for the selected features and filters candidates to just those targest that are within the region.

## Changes
 1. Remove code for dereplicating candidates in the in the LSR
 2. add code to set the candidate flag to false for any guide that has a hamming distance of  <=1 int the LSR  to any target sequence LSR in the full set (317 million for humans and SpCas9)
  - a  hardware popcount (.count_ones()) approach may be appropriate for rapidly screening 
 3. During interval screening a column should be added to the Polars target dataframe `feature_keys` with the foreign keys to the features dataframe
 4. make the default LSR 20
 5. I anticipate that it will be faster to do interval screening first to reduce the number of candidates that  hamming <=1 testing distance needs to be calulated for

## Data types

The input guide data will look like this after i add a primary_key column 

```
   candidate                   seq         chrom  start   stop  strand
       True   1567635393809842176  NC_000001.11  10445  10470    True
       True   6270548258482962432  NC_000001.11  10458  10483    True
       True  12418475214316847104  NC_000001.11  10471  10496    True
       True  12780412709848301568  NC_000001.11  10472  10497    True
       True  11918043498572906496  NC_000001.11  10484  10509    True
```

The input feature data will look like this 

```
primary_key      chrom  feature_start  feature_end  strand                 feature_id feature_type
0               NC_000001.11              0    248956422    True  NC_000001.11:1..248956422       region
1               NC_000001.11          10953        11523    True          gene-LOC102723747   pseudogene
2               NC_000001.11          10953        11523    True          id-LOC102723747-1         exon
3               NC_000001.11          11873        14409    True               gene-DDX11L1   pseudogene
4               NC_000001.11          11873        14409    True            rna-NR_046018.2   transcript
```