import os
import tempfile
import subprocess
import pytest
import pandas as pd
import numpy as np
from predictPAM import main

cwd = os.getcwd()
tempdir = os.path.join(cwd,"test_data")
gnbank = [os.path.join(tempdir, "Pseudomonas_aeruginosa_PAO1_107.gbk")]
pamseq='AGG'
strand = 'reverse'
threads = 1
targetlength = 25


main.get_fastas(gnbank, tempdir)

map_out = main.map_pam(tempdir, pamseq,threads,strand)

target_out = main.get_target(tempdir, map_out, targetlength, strand)
dict(list(target_out.items())[0:5])

parse_target_out = main.parse_target(target_out,strand,seqlengthtopam=12)
dict(list(parse_target_out.items())[0:5])

filterparsetargetdict= main.filter_parse_target(parse_target_out, threads=1, levendistance=2)

# reformating filterparsetargetdict- to make compatible for pybed
filterparsetargetdict_pd = pd.DataFrame.from_dict(filterparsetargetdict, orient='index')
# remove index, which is key of dict
filterparsetargetdict_pd_dindx = filterparsetargetdict_pd.reset_index(drop=True)
# pybed takes tab separated file with no header, plus first three column has to be as above
filterparsetargetdict_pd_dindx_tab = filterparsetargetdict_pd_dindx.to_csv(index=False,sep='\t',header=False)

genebankfeatures = main.get_genbank_features(gnbank)


# reformating genebankfeatures- to make compatible for pybed
genebankfeatures_df = pd.DataFrame(genebankfeatures)
# enpty tab crates isses in runnign pybed, so replance NaN with NA, then make tab separated
genebankfeatures_df = genebankfeatures_df.replace(np.nan, "NA")
genebankfeatures_df_tab = genebankfeatures_df.to_csv(index=False,sep='\t',header=False)


# pybedtools

down, up = main.get_nearby_feature(filterparsetargetdict_pd_dindx_tab, genebankfeatures_df_tab)


### Adding column name to upstream and downstream output from pybed.
# get columns from two input file: target_mappingfile(filterparsetargetdict), featurefile(genebankfeatures_df)
targetfile_columns = list(list(filterparsetargetdict.values())[0].keys())
featurefile_columns = list(genebankfeatures_df.columns) # from genebank feature dataframe
joined_columns = targetfile_columns + featurefile_columns
joined_columns.append("distance") # on top of mapping upstream and downstream, we are recording distance in pybed closest. Thus need an extra column


# merge upstream and downstream outputfilename
outfilename = os.path.join(tempdir,"outfile.txt")
merged_down_ups = main.merge_downstream_upstream(down,up,joined_columns,outputfilename=outfilename)
merged_down_ups.head()


merged_down_ups.to_csv('reverse_AGG.csv')
