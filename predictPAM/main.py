#!/usr/bin/env python
'''

'''
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
    parser.add_argument('--pamseq', '-p', type=str, required=True, help='A short PAM motif to search for, may be use IUPAC ambiguous alphabet')
    parser.add_argument('--outfile', '-o', type=str, required=True, help='The table of pam sites and data')
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    parser.add_argument('--keeptemp' ,help="Should intermediate files be kept?", action='store_true')
    parser.add_argument('--log' ,help="Log file", default="predictPAM.log")
    parser.add_argument('--threads' ,help="Number of processor threads to use.", type=int, default=1)
    return parser
    

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




############### MAPPING ###########################################
def map_pam(tempdir, pamseq, threads=1):
    """Runs seqkit locate to find the pam in given genome (FASTA)
    
    Args:
        threads (int or str):the number of processor threads to use
        
    """
    try:
        infasta = os.path.join(tempdir, "out.fasta")
        parameters = ["seqkit",
                       "locate",
                       "-p", pamseq,
                       infasta,"-P",
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
        
        
######## Retrieving target sequence ################################
def get_target(tempdir,mappingfile):
    mapping_list = []
    infasta = SeqIO.read(os.path.join(tempdir, "out.fasta"), "fasta")
    bylines = mappingfile.splitlines()
    for entry in bylines:
        tline = entry.split()
        start_ps = int(tline[1])
        pam_start = start_ps-1
        cut_ps = start_ps - 5
        target_sp = start_ps - 22
        target_ep = start_ps-2
        target_seq = str(infasta.seq)[target_sp:target_ep]
        strand = "+"
        accession = infasta.id
        pamseq = tline[3]
        #if PAM match at the beginning or the end then  target seq might be smaller or none
        if len(target_seq) == 20:
            mapping_list.append({"accession": accession, "target_sp": target_sp,
            "target_ep": target_ep, "target_seq": target_seq, "pam_seq": pamseq,
            "cut-pos": cut_ps, "pam-pos": pam_start, "strand": "+"})
        return mapping_list

#############
def get_cds(genebank_record):
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

    # map = BedTool("mapping.txt")
# map.head()
#
# # downstream of target sequence
# downstream = map.closest('features.txt', d=True, fd=True, D="a", t="first")
# downstream.head()
#
# # upstream of target sequence
# upstream = map.closest('features.txt', d=True, id=True, D="a", t="first")
# upstream.head()
#
# # add row names to the dataframe
# n = downstream.to_dataframe().shape[1]
# import string
# rownames = list(string.ascii_uppercase[0:n])
#
# # convert bed file to dataframe
# downstream_df = downstream.to_dataframe(names=rownames)
# downstream_df.to_csv('downstream.txt', sep='\t', index=False, header=True)
#
# upstream_df = upstream.to_dataframe(names=rownames)
# upstream_df.to_csv('upstream.txt', sep='\t', index=False, header=True)
#
# # merge downstream and upstream information and export as csv
# all_df = pd.merge(downstream_df, upstream_df,
#                   right_on=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
#                   left_on=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
#                   how='outer')
#
# all_df.to_csv('all.txt', sep='\t', index=False, header=True)


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
        SeqIO.write(SeqIO.parse(args.gbkfile, "genbank"), os.path.join(tempdir, "out.fasta"), "fasta")
        print(tempdir)
        # mapping
        logging.info("Mapping pam to the genome")
        mapfile = map_pam(tempdir=tempdir, pamseq=args.pamseq, threads=args.threads)
        
        #Retrieving target sequence
        logging.info("Retrieving target sequence for matching PAM : %s", tempdir)
        targetlist = get_target(tempdir=tempdir, mappingfile=mapfile)
        
        # get cds list
        logging.info("Retrieving CDS information")
        cdslist = get_cds(SeqIO.parse(args.gbkfile, "genbank"))
        
        # get gene list
        logging.info("Retrieving Gene information")
        genelist = get_gene(SeqIO.parse(args.gbkfile, "genbank"))
        
        # merge cds and genelist
        logging.info("Merge cda and gene file : %s", tempdir)
        merge_cds_gene(cdslist=cdslist,genelist=genelist, tempdir=tempdir)
        
        # pybedtools
        logging.info("Mapping features downstream of target sequence")
        find_downstream = pybed_downstream(tempdir,mapfile_from_pam=mapfile)
        
        logging.info("Mapping features upstream of target sequence")
        find_upstream = pybed_upstream(tempdir,mapfile_from_pam=mapfile)
        
        # merge upstream and downstream output
        logging.info("Writing features file : %s", tempdir)
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
