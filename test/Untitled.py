target_dict = {}
# track keys so that any duplicated entry can be removed from the final dictionary
keys_list =[]
# this won't work for big genomes because it reads into memory try seqio index
infasta = Fasta(os.path.join(tempdir, "forward.fasta"))
infasta_complement = Fasta(os.path.join(tempdir, "reverse.fasta"))
bylines = map_out.splitlines()
for entry in bylines:
    tline = entry.split()
    whichchromose=tline[0]
    pam_sp = int(tline[1]) - 1 # -1 to adjust- because seqkit convention - starts from 1 but in python starts from 0.
    pam_ep = int(tline[2])
    pam_seq = tline[3]
    seqid = infasta[whichchromose].name
    fastalength = infasta[whichchromose].unpadded_len
    # note the strand is not mean + from seqkit mapping. Comes from user- orientation of genome to search for target
    if strand=="forward":
        target_sp = pam_sp - targetlength ## this might give use negative value.
        target_ep = pam_sp
        if (target_sp > 0) and (target_ep > target_sp):
            target_seq = infasta[whichchromose][target_sp:target_ep].seq
    if strand=="reverse":
        target_sp = pam_ep
        target_ep = pam_ep + targetlength
        if (target_ep > target_sp) and (target_ep < fastalength):
            target_seq = infasta_complement[whichchromose][target_sp:target_ep].seq
    target_dict[target_seq]= {"seqid": seqid, "target_sp": target_sp,
            "target_ep": target_ep, "pam_seq": pam_seq,
             "strand": strand}
    keys_list.append(target_seq)

remove_target = [k for k, v in Counter(keys_list).items() if v > 1] # list of keys with more than one observation
target_dict2 = {k: v for k, v in target_dict.items() if k not in remove_target} # although dict over writes on non-unique key, but we want to complete remove such observation
return target_dict2


target_dict = {}
# track keys so that any duplicated entry can be removed from the final dictionary
keys_list =[]
# this won't work for big genomes because it reads into memory try seqio index
infasta = Fasta(os.path.join(tempdir, "forward.fasta"))
infasta_complement = Fasta(os.path.join(tempdir, "reverse.fasta"))
bylines = map_out.splitlines()
for entry in bylines:
    tline = entry.split()
    whichchromose=tline[0]
    pam_sp = int(tline[1]) - 1 # -1 to adjust- because seqkit convention - starts from 1 but in python starts from 0.
    pam_ep = int(tline[2])
    pam_seq = tline[3]
    seqid = infasta[whichchromose].name
    # note the strand is not mean + from seqkit mapping. Comes from user- orientation of genome to search for target
    if strand=="forward":
        target_sp = pam_sp - targetlength ## this might give use negative value.
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



filelist = [os.path.join(tempdir, "Pseudomonas_aeruginosa_PAO1_107.gbk")]

f1 = open(os.path.join(tempdir,"forward.fasta"),"w")
f2 = open(os.path.join(tempdir,"reverse_complement.fasta"),"w")
for file in filelist:
    recs = SeqIO.parse(file, "genbank")
    for rec in recs:
        record_f = rec
        record_rc = rec.reverse_complement()
        record_rc.id =rec.id
        record_rc.name = rec.name+"_reverse_complement"
        record_rc.description = rec.description+"_reverse_complement"
        SeqIO.write(record_f,f1,"fasta")
        SeqIO.write(record_rc,f2,"fasta")



def get_fastas(filelist, tempdir):
    """Returns Fasta and complement of Fasta for a given Genbank file

    Args:
        filelist (str): Genbank file to process

    Returns:
        forward.fasta(str): Fasta file in forward orientation (5'-3')
        reverse_complement.fasta(str): Reverse Complement of Fasta file
    """
    try:
        f1 = open(os.path.join(tempdir,"forward.fasta"),"w")
        f2 = open(os.path.join(tempdir,"reverse_complement.fasta"),"w")
        for file in filelist:
            recs = SeqIO.parse(file, "genbank")
            for rec in recs:
                record_f = rec
                record_rc = rec.reverse_complement()
                record_rc.id =rec.id
                record_rc.name = rec.name+"_reverse_complement"
                record_rc.description = rec.description+"_reverse_complement"
                SeqIO.write(record_f,f1,"fasta")
                SeqIO.write(record_rc,f2,"fasta")
    except Exception as e:
        print("An error occurred in input genbank file")
        raise e



filelist = [os.path.join(tempdir, "Pseudomonas_aeruginosa_PAO1_107.gbk")]
