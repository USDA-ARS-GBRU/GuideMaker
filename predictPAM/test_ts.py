genbank = "Burkholderia_thailandensis_E264__ATCC_700388_133.gbk"
pamseq ="ATCGA"
tempdir = '/var/folders/52/rbrrfj5d369c35kd2xrktf3m0000gq/T/pamPredict_cs6mazl9'
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

filterparsetargetdict = filter_parse_out
# reformating filterparsetargetdict- to make compatible for pybed
filterparsetargetdict_pd = pd.DataFrame.from_dict(filterparsetargetdict, orient='index')
# remove index, which is key of dict
filterparsetargetdict_pd_dindx = filterparsetargetdict_pd.reset_index(drop=True)
# pybed takes tab separated file with no header, plus first three column has to be as above
filterparsetargetdict_pd_dindx_tab = filterparsetargetdict_pd_dindx.to_csv(index=False,sep='\t',header=False)


############
genebankfeatures = get_genbank_features(genbank)
# reformating genebankfeatures- to make compatible for pybed
genebankfeatures_df = pd.DataFrame(genebankfeatures)
# enpty tab crates isses in runnign pybed, so replance NaN with NA, then make tab separated
genebankfeatures_df = genebankfeatures_df.replace(np.nan, "NA")
genebankfeatures_df_tab = genebankfeatures_df.to_csv(index=False,sep='\t',header=False)

down, up = get_nearby_feature(filterparsetargetdict_pd_dindx_tab, genebankfeatures_df_tab)

### add column name to
kk = filterpasrsedict.values()
mm=list(kk)
nn=mm[0]
oo= list(nn.keys())
pp = list(gnbfeature_df.columns)
joined = oo+pp
joined.append("distance")


mm = merge_downstream_upstream(down,up,joined)



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




merge_downstream_upstream(downsfile,upsfile,columns_name)


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


targetfile_columns = list(list(filterparsetargetdict.values())[0].keys())
featurefile_columns = list(gnbfeature_df.columns)
joined_columns = targetfile_columns + featurefile_columns
joined_columns.append("distance")
#######################################################
##############3 parding siqio.index ####################
genbank = "Burkholderia_thailandensis_E264__ATCC_700388_133.gbk"
record_index = SeqIO.index(genbank, "genbank")
ids = list(record_index)
records_forward = (record_index[fid] for fid in ids)
records_reverse = (SeqRecord(seq=record_index[rid].seq.complement(),
                            id=record_index[rid].id,
                            description=record_index[rid].description+"_complement",
                            name=record_index[rid].name+"_complement"
                            ) for rid in ids)
SeqIO.write(records_forward, "out.fasta", "fasta")
SeqIO.write(records_reverse, "out_complement.fasta", "fasta")


# get list of recode in gbkfile
list(gbk_indx.keys())

#
list(gbk_indx.values())

# get values for first record
v = list(gbk_indx.values())[0]
type(v)

# get seqeuce
v.seq


tempdir='/var/folders/52/rbrrfj5d369c35kd2xrktf3m0000gq/T/pamPredict_z8dgknrz'


ids = list(record_index)
a=ids[0]
records = (record_index[id] for id in ids)

ids = list(record_index)
records_reverse = (SeqRecord(seq=record_index[rid].seq.complement(),
                            id=record_index[rid].id,
                            description=record_index[rid].description,
                            name=record_index[rid].name
                            ) for rid in ids)
SeqIO.write(records_reverse, "out_complement.fasta", "fasta")


tempdir='/var/folders/71/l663_yxn40n1f2l111k9x7nc0000gn/T/pamPredict__k9r309z'
def map_pam(tempdir, pamseq, threads, strand):
    """Runs seqkit locate to find the PAM site in the genome (Fasta)

    Args:
        threads (int): the number of processor threads to use

    Returns:
        (str): Data in Bedfile format with matches

    """

infasta = os.path.join(tempdir, "out.fasta")
infasta_complement = os.path.join(tempdir, "out_complement.fasta")
parameters = ["seqkit",
               "locate",
               "-p", pamseq,
               infasta,"-P",
               "--bed",
                "--threads", str(threads),
                "-d"]
if strand == "reverse":
    parameters[4] =infasta_complement
p = subprocess.Popen(parameters, stdout=subprocess.PIPE)
out, err = p.communicate()
out = out.decode('utf-8')

mappingdata=out
target_dict = {}
# track keys so that any duplicated entry can be removed from the final dictionary
keys_list =[]
# this won't work for big genomes because it reads into memory try seqio index
infasta = SeqIO.index(os.path.join(tempdir, "out.fasta"), "fasta")
infasta_complement = SeqIO.index(os.path.join(tempdir, "out_complement.fasta"), "fasta")
bylines = mappingdata.splitlines()
for entry in bylines:
    tline = entry.split()
    whichchromose=tline[0]
    pam_sp = int(tline[1]) - 1 # -1 to adjust- because seqkit convention - starts from 1 but in python starts from 0.
    pam_ep = int(tline[2])
    pam_seq = tline[3]
    seqid = infasta[whichchromose].id
    # note the strand is not mean + from seqkit mapping. Comes from user- orientation of genome to search for target
    if strand=="forward":
        target_sp = pam_sp - targetlength
        target_ep = pam_sp
        target_seq = str(infasta[whichchromose].seq)[target_sp:target_ep]
    if strand=="reverse":
        target_sp = pam_ep
        target_ep = pam_ep + targetlength
        target_seq = str(infasta_complement[whichchromose].seq)[target_sp:target_ep]
    if len(target_seq) == targetlength:
            target_dict[target_seq]= {"seqid": seqid, "target_sp": target_sp,
            "target_ep": target_ep, "pam_seq": pam_seq,
             "strand": strand}
    keys_list.append(target_seq)

remove_target = [k for k, v in Counter(keys_list).items() if v > 1] # list of keys with more than one observation
target_dict2 = {k: v for k, v in target_dict.items() if k not in remove_target} # although dict over writes on non-unique key, but we want to complete remove such observation


aa = get_target(tempdir, out, targetlength, strand)


##########3 learnign pyfaidx

from pyfaidx import Fasta
genes = Fasta(os.path.join(tempdir, "out.fasta"))
genes.keys()
genes['NC_007650'][200:300].seq


mappingdata=map_out
target_dict = {}
# track keys so that any duplicated entry can be removed from the final dictionary
keys_list =[]
# this won't work for big genomes because it reads into memory try seqio index
infasta = Fasta(os.path.join(tempdir, "out.fasta"))
infasta_complement = Fasta(os.path.join(tempdir, "out_complement.fasta"))
bylines = mappingdata.splitlines()
for entry in bylines:
    tline = entry.split()
    whichchromose=tline[0]
    pam_sp = int(tline[1]) - 1 # -1 to adjust- because seqkit convention - starts from 1 but in python starts from 0.
    pam_ep = int(tline[2])
    pam_seq = tline[3]
    seqid = infasta[whichchromose].name
    # note the strand is not mean + from seqkit mapping. Comes from user- orientation of genome to search for target
    if strand=="forward":
        target_sp = pam_sp - targetlength
        target_ep = pam_sp
        target_seq = infasta[whichchromose][target_sp:target_ep].seq
    if strand=="reverse":
        target_sp = pam_ep
        target_ep = pam_ep + targetlength
        target_seq = infasta_complement[whichchromose][target_sp:target_ep].seq
    if len(target_seq) == targetlength:
            target_dict[target_seq]= {"seqid": seqid, "target_sp": target_sp,
            "target_ep": target_ep, "pam_seq": pam_seq,
             "strand": strand}
    keys_list.append(target_seq)

remove_target = [k for k, v in Counter(keys_list).items() if v > 1] # list of keys with more than one observation
target_dict2 = {k: v for k, v in target_dict.items() if k not in remove_target} # although dict over writes on non-unique key, but we want to complete remove such observation












#########

def get_fastas(*args):
    f_list=[]
    r_list=[]
    for file in args:
        record_index = SeqIO.index(file, "genbank")
        ids = list(record_index)
        records_forward = [record_index[fid] for fid in ids]
        f_list.append(records_forward)
        records_reverse = [SeqRecord(seq=record_index[rid].seq.complement(),
                                    id=record_index[rid].id,
                                    description=record_index[rid].description+"_complement",
                                    name=record_index[rid].name+"_complement"
                                    ) for rid in ids]
        r_list.append(records_reverse)
    return f_list, r_list
                                    

        SeqIO.write(records_forward , "out.fasta", "fasta")
        SeqIO.write(records_reverse, "out_complement.fasta", "fasta")


aa =get_fastas("Pseudomonas_aeruginosa_PAO1_107.gbk","Burkholderia_thailandensis_E264__ATCC_700388_133.gbk")



def get_fastas(*args):
    f_list=[]
    for file in args:
        record_index = SeqIO.index(file, "genbank")
        ids = list(record_index)
        records_forward = [record_index[fid] for fid in ids]
        f_list.append(records_forward)
    return f_list


aa = get_fastas("Pseudomonas_aeruginosa_PAO1_107.gbk","Burkholderia_thailandensis_E264__ATCC_700388_133.gbk")
SeqIO.write(aa , "out.fasta", "fasta")


record_index = SeqIO.index("Burkholderia_thailandensis_E264__ATCC_700388_133.gbk", "genbank")
record_index_list=list(record_index)
f_list=[]
for item in record_index_list:
    record= SeqRecord(record_index[item].seq,record_index[item].id,
            record_index[item].name,
            record_index[item].description)
    f_list.append(record)

SeqIO.write(f_list , "out.fasta", "fasta")
    
    
    
records = (SeqRecord(Seq(seq, generic_dna), str(index)) for index,seq in enumerate(sequence_set) )
SeqIO.write(records, file_location, "fasta")



##
from Bio import SeqIO
for index, record in enumerate(SeqIO.parse("Burkholderia_thailandensis_E264__ATCC_700388_133.gbk", "genbank")):
    print("index %i, ID = %s, length %i, with %i features"
          % (index, record.id, len(record.seq), len(record.features)))


##
def get_fastas(*args):
    f_list=[]
    r_list=[]
    for file in args:
        record_index = SeqIO.index(file, "genbank")
        record_index_list=list(record_index)
        for item in record_index_list:
            record_f = SeqRecord(record_index[item].seq,
                              record_index[item].id,
                              record_index[item].name,
                              record_index[item].description)
            record_r = SeqRecord(record_index[item].seq.complement(),
                              record_index[item].id+"_complement",
                              record_index[item].name+"_complement",
                              record_index[item].description+"_complement")
            f_list.append(record_f)
            r_list.append(record_r)
    SeqIO.write(f_list , "out.fasta", "fasta")
    SeqIO.write(r_list , "out_complement.fasta", "fasta")


# def get_fastas(*args):
#     f_list=[]
#     r_list=[]
#     for file in args:
#         record_index = SeqIO.index(file, "genbank")
#         record_index_list=list(record_index)
#         for item in record_index_list:
#             record_f = SeqRecord(record_index[item].seq,
#                               record_index[item].id)
#             record_r = SeqRecord(record_index[item].seq.complement(),
#                               record_index[item].id)
#             f_list.append(record_f)
#             r_list.append(record_r)
#     SeqIO.write(f_list , "out.fasta", "fasta")
#     SeqIO.write(r_list , "out_complement.fasta", "fasta")



get_fastas("Pseudomonas_aeruginosa_PAO1_107.gbk","Burkholderia_thailandensis_E264__ATCC_700388_133.gbk")



def get_fastas(*args):
    f1 = open("example.fasta", "w")
    for file in args:
        records = SeqIO.parse(file, "genbank")
        SeqIO.write(records, f1, "fasta")
        

############
file="Burkholderia_thailandensis_E264__ATCC_700388_133.gbk"
record = SeqIO.parse(file, "genbank")
SeqIO.write(record, "test.out", "fasta")
  


recs = SeqIO.parse(sys.argv[1], 'fastq')

for rec in recs:
    rc = rec.reverse_complement()
    rc.id = rec.id
    rc.name = rec.name + '_R2.fastq'
    rc.description = ''
    SeqIO.write(rc, rc.name, 'fastq')
    
def get_fastas(*args):
    f1 = open("forward.fasta", "w")
    f2 = open("reverse.fasta", "w")
    for file in args:
        recs = SeqIO.parse(file, "genbank")
        for rec in recs:
            record_f = rec
            record_r = SeqRecord(rec.seq.complement(),rec.id,rec.name+"_complement",rec.description+"_complement")
            SeqIO.write(record_f, f1, "fasta")
            SeqIO.write(record_r,f2, "fasta")
        
get_fastas("Pseudomonas_aeruginosa_PAO1_107.gbk","Burkholderia_thailandensis_E264__ATCC_700388_133.gbk")



def get_fastas(*args):
    f1 = open(os.path.join(tempdir,"forward.fasta"),"w")
    f2 = open(os.path.join(tempdir,"reverse.fasta"),"w")
    for file in args:
        recs = SeqIO.parse(file, "genbank")
        for rec in recs:
            record_f = rec
            record_r = SeqRecord(rec.seq.complement(),rec.id,rec.name+"_complement",rec.description+"_complement")
            SeqIO.write(record_f , f1,"fasta")
            SeqIO.write(record_r , f2, "fasta")
