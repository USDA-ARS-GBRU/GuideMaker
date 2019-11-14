genbank = "Pseudomonas_aeruginosa_PAO1_107.gbk"
pamseq ="ATCGA"
tempdir ='/var/folders/71/l663_yxn40n1f2l111k9x7nc0000gn/T/pamPredict_91ot3oyt'
threads=2
strand = "forward"
targetlength = 25
seqlengthtopam = 12
levendistance = 2
get_fastas(genbank, tempdir)
map_out = map_pam(tempdir, pamseq, threads, strand)
target_out = get_target(tempdir, map_out, targetlength, strand)
list(target_out.items())[:4]

parse_out = parse_target(target_out, strand, seqlengthtopam)
list(parse_out.items())[:4]


filter_parse_out = filter_parse_target(parse_out, threads, levendistance)

list(filter_parse_out.items())[:4]


reformat_map = reformat_parse_target_for_pybed(filter_parse_out)

mapfile_from_pam = reformat_map
mapbed = BedTool(mapfile_from_pam.splitlines())


gnbfeature = get_genbank_features(genbank)
gnbfeature_df = pd.DataFrame(gnbfeature)
# enpty tab crates isses in runnign pybed, so replance NaN with NA
gnbfeature_df = gnbfeature_df.replace(np.nan, "NA")
gnbfeature_df.to_csv("test.csv")
gnbfeature_df_tab = gnbfeature_df.to_csv(index=False,sep='\t',header=False)
gnbfeature_df_tab = gnbfeature_df_tab.splitlines()
featurebed = BedTool(gnbfeature_df_tab)
downstream = mapbed.closest(featurebed , d=True, fd=True, D="a", t="first")
upstream = mapbed.closest(featurebed , d=True, id=True, D="a", t="first")

### add column name to
kk = filterpasrsedict.values()
mm=list(kk)
nn=mm[0]
oo= list(nn.keys())
pp = list(gnbfeature_df.columns)
joined = oo+pp
joined.append("distance")

## pybed to dataframe
downstream_df = downstream.to_dataframe(names=joined,low_memory=False)
downstream_df.to_csv("downstream_df.csv")
upstream_df = upstream.to_dataframe(names=joined,low_memory=False)
upstream_df.to_csv("upstream_df.csv")
all_df = pd.merge(downstream_df, upstream_df,
                  right_on=joined[:8],
                  left_on=joined[:8],
                  how='outer')
all_df.to_csv("all_df.csv")



### merging ups and downsfile
n = downsfile.to_dataframe().shape[1]
rownames = list(string.ascii_uppercase[0:n])
downstream_df = downsfile.to_dataframe(names=rownames,low_memory=False)
upstream_df = upsfile.to_dataframe(names=rownames,low_memory=False)
all_df = pd.merge(downstream_df, upstream_df,
                  right_on=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                  left_on=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                  how='outer')
all_df.to_csv(os.path.join(tempdir, "all.txt"), sep='\t', header=True, index=False)







dict={}
for l in featurebed:
    bedl = BedTool(l)
    dict["down"] = mapbed.closest(bedl, d=True, fd=True, D="a", t="first")
    dict["up"] = mapbed.closest(bedl , d=True, id=True, D="a", t="first")






downstream = featurebed.closest(l , d=True, fd=True, D="a", t="first")


gnbfeature_df_tab = gnbfeature_df.to_csv(index=False,sep='\t',header=False)
aa = gnbfeature_df_tab.splitlines()
test_list={}
for entry in aa:
    featurebed = BedTool(entry)
    print(featurebed)



downstream_df = downsfile.to_dataframe(low_memory=False)

coln = list(gnbfeature_df.columns)
downstream_df.shape











feature_list = []
genebank_file = SeqIO.parse(genbank,"genbank")
for entry in genebank_file:
    for record in entry.features:
        feature_dict = {}
        if record.type in ['CDS', 'gene']:
            feature_dict["accession"] = entry.id
            feature_dict["start"] = record.location.start.position
            feature_dict["stop"] = record.location.end.position
            feature_dict["type"] = record.type
            feature_dict["strand"] = 'reverse' if record.strand < 0 else 'forward'
            for qualifier_key, qualifier_val in record.qualifiers.items():
                feature_dict[qualifier_key] = qualifier_val
            feature_list.append(feature_dict)

len(feature_list)
