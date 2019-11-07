#!/usr/bin/env python
"""
predictPAM: A python module to predict custom PAM sites in any small genome

"""
# Author: Ravin Poudel

import sys
import os
import logging
import shutil
import random
import itertools
import subprocess
import matplotlib
import numpy
import tempfile
import argparse
import string
import pickle
import nmslib
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from itertools import permutations
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pybedtools import BedTool

def myparser():
    parser = argparse.ArgumentParser(description='predictPAM: A python module to predict custom PAM sites in any small genome')
    parser.add_argument('--gbkfile', '-i', type=str, required=True,
                        help='A genbank .gbk file')
    parser.add_argument('--pamseq', '-p', type=str, required=True, help='A short PAM motif to search for, may be use IUPAC ambiguous alphabet'),
    parser.add_argument('--targetlength', '-l', type=int, default=22, help='Length of the target sequence'),
    parser.add_argument('--strand', '-s', type=str, default="+", help='Strand of DNA'),
    parser.add_argument('--outfile', '-o', type=str, required=True, help='The table of pam sites and data')
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    parser.add_argument('--keeptemp' ,help="Should intermediate files be kept?", action='store_true')
    parser.add_argument('--log' ,help="Log file", default="predictPAM.log")
    parser.add_argument('--threads' ,help="Number of processor threads to use.", type=int, default=1)
    return parser

# GLOBAL parsing variables
ALLOWED_FILE_EXTENSIONS = set(['gbk'])


def _logger_setup(logfile):
    """Set up logging to a logfile and the terminal standard out.

    Args:
        fastq (str): The path to a fastq or fastq.gz file
        fastq2 (str): The path to a fastq or fastq.gz file for the reverse sequences

    """
    try:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%m-%d %H:%M',
                            filename=logfile,
                            filemode='w')
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        logging.getLogger('').addHandler(console)
    except Exception as e:
        print("An error occurred setting up logging")
        raise e

####################
def allowed_file(filename):
    """ Check for the correct file type
    
    Args:
        filename(file): File to process
    
    Returns:
        Booleans(): True or False
    
    """
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_FILE_EXTENSIONS
    
    

def get_fastas(genbank, tempdir):
    """Retruns fasta and complement of fasta for a given genbank file

    Args:
        genbank (genebank): Genbank file to process

    Returns:
        (out.fasta): Fasta file in forward orientation (5'-3')
        (out_complement.fasta): Complement of Fasta file in reverese orientation (3'-5')

    """
    if allowed_file(genbank):
        sequences = []
        sequences_complement=[]
        for record in SeqIO.parse(genbank, "genbank"):
            sequences.append(record)
            sequences_complement.append(SeqRecord(record.seq.complement(), record.id+"_complement",
            description=record.description+"_complement", name=record.name+"_complement"))
        SeqIO.write(sequences, os.path.join(tempdir,"out.fasta"), "fasta")
        SeqIO.write(sequences_complement, os.path.join(tempdir,"out_complement.fasta"), "fasta")
    else:
        print("Invalid file type, allowed file type is .gbk")
        

############### MAPPING ###########################################
def map_pam(tempdir, pamseq,threads,strand):
    """Runs seqkit locate to find the pam in given genome (FASTA)

    Args:
        threads (int or str):the number of processor threads to use

    Returns:
        (str): Bedfile format with matches

    """
    try:
        infasta = os.path.join(tempdir, "out.fasta")
        infasta_complement = os.path.join(tempdir, "out_complement.fasta")
        if strand == "+":
            parameters = ["seqkit",
                           "locate",
                           "-p", pamseq,
                           infasta,"-P",
                           "--bed",
                            "--threads", str(threads)]
            p = subprocess.Popen(parameters, stdout=subprocess.PIPE)
            out, err = p.communicate()
            out = out.decode('utf-8')
        if strand =="-":
            parameters = ["seqkit",
                           "locate",
                           "-p", pamseq,
                           infasta_complement,"-P",
                           "--bed",
                            "--threads", str(threads)]
            p = subprocess.Popen(parameters, stdout=subprocess.PIPE)
            out, err = p.communicate()
            out = out.decode('utf-8')
        return out
    except subprocess.CalledProcessError as e:
        raise e
    except FileNotFoundError as f:
        raise f

tempdir='/var/folders/52/rbrrfj5d369c35kd2xrktf3m0000gq/T/pamPredict_9z8641gp'
######## Retrieving target sequence ################################
def get_target(tempdir, mappingdata, targetlength, strand):
    """Given bedfile of PAM locations, function goes up or down and finds target

    Args:
        tempdir (str): Temporary directory
        mappingdata (str): Bedfile output from map_pam

    Returns:
        (str): Bedfile format with matches
    
    """
    target_dict = {}
    infasta = SeqIO.read(os.path.join(tempdir, "out.fasta"), "fasta")
    infasta_complement = SeqIO.read(os.path.join(tempdir, "out_complement.fasta"), "fasta")
    bylines = mappingdata.splitlines()
    if strand == "+":
        for entry in bylines:
            tline = entry.split()
            pam_start_ps = int(tline[1]) - 1 # -1 to adjust- beacaue seqkit convention - starts from 1 but in python starts from 0.
            pam_end_ps = int(tline[2])
            pam_seq = tline[3]
            target_sp = pam_start_ps - targetlength
            target_seq = str(infasta.seq)[target_sp:pam_start_ps]
            seqid = infasta.id
            #if PAM match at the beginning or the end then  target seq might be smaller or none
            if len(target_seq) == targetlength:
                target_dict[target_seq]= {"seqid": seqid, "target_sp": target_sp,
                "target_ep": pam_start_ps, "pam_seq": pam_seq,
                "pam_start_ps": pam_start_ps, "strand": strand}
    if strand == "-":
        for entry in bylines:
            tline = entry.split()
            pam_start_ps = int(tline[1]) - 1 # -1 to adjust- beacaue seqkit convention - starts from 1 but in python starts from 0.
            pam_end_ps = int(tline[2])
            pam_seq = tline[3]
            target_ep = pam_end_ps + targetlength
            target_seq = str(infasta_complement.seq)[pam_end_ps:target_ep]
            seqid = infasta.id
            #if PAM match at the beginning or the end then  target seq might be smaller or none
            if len(target_seq) == targetlength:
                target_dict[target_seq]={"seqid": seqid, "target_sp": pam_end_ps,
                "target_ep": target_ep, "pam_seq": pam_seq,
                "pam_start_ps": pam_start_ps, "strand": strand}
    return target_dict
############### Means to filter target sequences ############
# parse key then create 12 bp and remanining part, also take care of strand specifictity
def parse_target(targetdict,strand):
    """Given a dictionary of target sequence, parse target sequence into two parts:
    close12 and remainingseq, then create a new dictionary with unique close12 sequences as keys

        Args:
            targetdict (dict): Dictionary of target sequences obtained after running get_target

        Returns:
            parse_target_dict(dict): Dictionary with unique close12 sequences as keys
        
    """
    parse_target_dict={}
    for items in targetdict.items():
        if items[1]['strand']=="+":
            close12 = items[0][-12:]
            remainingseq = items[0][:-12]
            items[1]['target'] = items[0] # move target sequence as value
            items[1]['remainingseq'] = remainingseq
            parse_target_dict[close12] = items[1] ## will retain only target with unique close12 bp
        if strand == "-":
            close12 = target_key[:12]
            remainingseq = target_key[12:]
            items[1]['target'] = items[0]
            items[1]['remainingseq'] = remainingseq
            parse_target_dict[close12] = items[1]
    return parse_target_dict

#nms lib - create a index with remainingseq
def create_index(strings):
    """Initializes and returns a NMSLIB index
    Args:
        strings (str): Strings to calculate distance
        
    Returns:
        index(index): Returns a NMSLIB index
    """
    index = nmslib.init(space='leven',
                        dtype=nmslib.DistType.INT,
                        data_type=nmslib.DataType.OBJECT_AS_STRING,
                        method='small_world_rand')
    index.addDataPointBatch(strings)
    index.createIndex(print_progress=True)
    return index

## Calculate Levenshtein distance among remainingseq, and remove any remainingseq that are similar
def filter_parse_target(parse_dict):
    """Returns a filter target sequences based on Leven distance (greater than 2) on sequences 12 bp away from the PAM motif
    Args:
        parse_dict(dict): A dictionary with parse target sequence
        
    Returns:
        filter_pasrse_dict(dict): A dictionary whose target sequences that are 12 bp away from PAM sequences are at the Levenshtein distance greater than 2.
    """
    filter_parse_dict={}
    # initialize a new index
    ref_index = create_index(list(parse_dict.keys()))
    for keys, value in parse_dict.items():
        ids, distances = ref_index.knnQuery(keys, k=100) ## k =number of k nearest neighbours (knn)
        check_values=[0,1,2] ## Levenshtein Distance greater than 2 is selected, less than 2 is considered as too similar -- maximizing specificity of target
        check = np.isin(distances, value) # if any element of value in the array, then True else False
        sum_check = np.sum(check) # If non present, all are false, hence sum == zero, else sum >=1
        if sum_check == 0 :
            filter_parse_dict[keys] = value
    return filter_parse_dict

def reformat_parse_target_for_pybed(filterpasrsedict):
    """Converts dictionary to tab separated format as need for pybed tools.
    Args:
        filterpasrsedict(dict): A dictionary with filter parse target sequence
        
    Returns:
        tab(file): Tab separated file
    """
    dict_to_pd = pd.DataFrame.from_dict(filterpasrsedict,orient='index')
    dict_to_pd.reset_index(inplace=True) # reindex
    dict_to_pd_reorder=dict_to_pd[['seqid', 'target_sp', 'target_ep', 'pam_seq', 'pam_start_ps','strand', 'target', 'remainingseq','index']] ## pybed need first column to be seqid/chromosome, next column start and thired column as end position
    dict_to_pd_reorder.columns = ['seqid', 'target_sp', 'target_ep', 'pam_seq', 'pam_start_ps',
    'strand', 'target', 'remainingseq','close12'] # rename index to close12
    dict_to_pd_reorder_tab = dict_to_pd_reorder.to_csv(index=False,sep='\t',header=False)
    return dict_to_pd_reorder_tab


###################
# filter_parse_dict need to convert to list with appropriate order then can be passed for pybed tools

#############
def get_cds(genebank_record):
    """Return a list of CDS's for a genbank file

    Args:
        genbank_record (obj): Biopython genbank object


    Returns:
        (list): List of CDS
    """
    cds_list = []
    for cds_record in genebank_record:
        for cds_feature in cds_record.features:
            if cds_feature.type == "CDS":
                start = cds_feature.location.start.position
                stop = cds_feature.location.end.position
                accession = cds_record.id
                des = cds_record.description
                type = cds_feature.type
                try:
                    geneID = cds_feature.qualifiers['db_xref'][0]
                except KeyError:
                    geneID = None
                try:
                    locus_tag = cds_feature.qualifiers['locus_tag'][0]
                except KeyError:
                    locus_tag = None
                try:
                    amino = cds_feature.qualifiers['translation'][0]
                except KeyError:
                    amino = None
                try:
                    product = cds_feature.qualifiers['product'][0]
                except KeyError:
                    product = None
                try:
                    other = cds_feature.qualifiers['note'][0]
                except KeyError:
                    other = None
                if cds_feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                cds_list.append({"accession": accession,  "start": start, "stop": stop,
                                 "strand": strand, "locus_tag": locus_tag, "type": type,
                                 "geneID": geneID, "product": product})
    return cds_list



def get_gene(genebank_record):
    gene_list = []
    for gene_record in genebank_record:
        for gene_feature in gene_record.features:
            if gene_feature.type == 'gene':
                start = gene_feature.location.start.position
                stop = gene_feature.location.end.position
                accession = gene_record.id
                try:
                    des = gene_feature.description
                except AttributeError:
                    des = None
                type = gene_feature.type
                try:
                    geneID = gene_feature.qualifiers['db_xref'][0]
                except KeyError:
                    geneID = None
                try:
                    locus_tag = gene_feature.qualifiers['locus_tag'][0]
                except KeyError:
                    locus_tag = None
                try:
                    gene_name = gene_feature.qualifiers['name'][0]
                except KeyError:
                    gene_name = None
                if gene_feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                gene_list.append({"accession": accession, "start": start, "stop": stop,
                                  "strand": strand, "locus_tag": locus_tag, "type": type, "gene_name": gene_name,
                                  "geneID": geneID})
    return gene_list

def merge_cds_gene(cdslist, genelist, tempdir):
    cds_df = pd.DataFrame(cdslist)
    gene_df = pd.DataFrame(genelist)
    merge_df = pd.merge(gene_df, cds_df,  on=["locus_tag"], how='outer')
    merge_df.to_csv(os.path.join(tempdir, "features.txt"), sep='\t', header=False, index=False)

########################################################################################################
# ######### pybedtools ########

def pybed_downstream(tempdir, mapfile_from_pam):
    featurefile = os.path.join(tempdir, "features.txt")
    mapbed = BedTool(mapfile_from_pam.splitlines())
    downstream = mapbed .closest(featurefile , d=True, fd=True, D="a", t="first")
    return downstream


def pybed_upstream(tempdir,mapfile_from_pam):
    featurefile = os.path.join(tempdir, "features.txt")
    mapbed = BedTool(mapfile_from_pam.splitlines())
    upstream = mapbed.closest(featurefile , d=True, id=True, D="a", t="first")
    return upstream


def merge_downstream_upstream(downsfile,upsfile,tempdir):
    n = downsfile.to_dataframe().shape[1]
    rownames = list(string.ascii_uppercase[0:n])
    downstream_df = downsfile.to_dataframe(names=rownames,low_memory=False)
    upstream_df = upsfile.to_dataframe(names=rownames,low_memory=False)
    all_df = pd.merge(downstream_df, upstream_df,
                      right_on=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                      left_on=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                      how='outer')
    all_df.to_csv(os.path.join(tempdir, "all.txt"), sep='\t', header=True, index=False)

################
def main(args=None):
    """Run Complete predictPAM workflow.
    """
    # Set up logging
    parser = myparser()
    if not args:
        args = parser.parse_args()

    _logger_setup(args.log)
    try:
        pamseq = Seq(args.pamseq, IUPAC.unambiguous_dna)
        #parse genbank file: todo allow multiple Genbank files
        #record = SeqIO.parse(args.gbkfile, "genbank")

        if args.tempdir:
            if not os.path.exists(args.tempdir):
                logging.warning("Specified location for tempfile ({}) does not exist, using default location.".format(tempdir))
                tempdir = tempfile.mkdtemp(prefix='pamPredict_')
        else:
            tempdir = tempfile.mkdtemp(prefix='pamPredict_', dir=args.tempdir)
        logging.info("Temp directory is: %s", tempdir)
        # Try to avoid writing out and reading genome in fasta format
        logging.info("Retriving fastas- forward and reverse from a genbanke file")
        get_fastas(genbank=args.gbkfile, tempdir=tempdir)
        print(tempdir)
        # mapping
        logging.info("Mapping pam to the genome")
        mapfile = map_pam(tempdir=tempdir, pamseq=args.pamseq, threads=args.threads, strand=args.strand)

        #Retrieving target sequence
        logging.info("Retrieving target sequence for matching PAM : %s", tempdir)
        targetdict = get_target(tempdir=tempdir, mappingdata=mapfile, targetlength=args.targetlength, strand=args.strand)
        
        # Parsing target sequence into two: 1)close12 and 2) remainingseq
        logging.info("Parsing target sequence into two: 1)close12 and 2) remainingseq")
        parsetargetdict = parse_target(targetdict, strand=args.strand)
        
        # Calculate Levenshtein distance among remainingseq, and remove any remainingseq that are similar
        logging.info("Filtering parse target sequence based on Levenshtein distance using NMSLIB index, with knn=100")
        filterparsetargetdict = filter_parse_target(parsetargetdict)
        
        ## Formatiing filter parse target sequnece as needed for pybed tools
        logging.info("Formatiing filter parse target sequnece as needed for pybed tools")
        tabfile_for_pybed= reformat_parse_target_for_pybed(filterparsetargetdict)


        # get cds list
        logging.info("Retrieving CDS information")
        cdslist = get_cds(SeqIO.parse(args.gbkfile, "genbank"))
        
        # get gene list
        logging.info("Retrieving Gene information")
        genelist = get_gene(SeqIO.parse(args.gbkfile, "genbank"))
        
        # merge cds and genelist
        logging.info("Merge cda and gene file- features.txt located at : %s", tempdir)
        merge_cds_gene(cdslist=cdslist,genelist=genelist, tempdir=tempdir)
        
        # pybedtools
        logging.info("Mapping features downstream of target sequence")
        find_downstream = pybed_downstream(tempdir,mapfile_from_pam=tabfile_for_pybed)
        
        logging.info("Mapping features upstream of target sequence")
        find_upstream = pybed_upstream(tempdir,mapfile_from_pam=tabfile_for_pybed)
        
        # merge upstream and downstream output
        logging.info("Writing features file all.txt - : %s", tempdir)
        merge_downstream_upstream(find_downstream,find_upstream,tempdir)

    except Exception as e:
        logging.error("predictPAM terminated with errors. See the log file for details.")
        logging.error(e)
        raise SystemExit(1)
    finally:
        try:
            if not args.keeptemp:
                shutil.rmtree(tempdir)
        except UnboundLocalError:
            pass
        except AttributeError:
            pass



if __name__ == '__main__':
    main()
