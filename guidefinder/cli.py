"""GuideFinder: The command line interface

"""
import logging
import argparse
import tempfile
import shutil
import os
import pybedtools

from Bio import SeqIO

import guidefinder


def myparser():
    parser = argparse.ArgumentParser(description='GuideFinder: globally design gRNAs for any CRISPR-Cas system in any small genome')
    parser.add_argument('--gbkfile', '-i', nargs='+', type=str, required=True, help='One or more genbank .gbk files for a single genome')
    parser.add_argument('--pamseq', '-p', type=str, required=True, help='A short PAM motif to search for, may be use IUPAC ambiguous alphabet')
    parser.add_argument('--pam_orientation', '-r', choices=['5prime', '3prime'], default='5prime', help="PAM position relative to target: 5prime: [PAM][target], 3prime: [target][PAM], Cas9 for example is 5prime")
    parser.add_argument('--targetlength', '-l', type=int, default=22, help='Length of the target sequence')
    parser.add_argument('--strand', '-s', choices=['forward','reverse', 'both'], default='both', help='Strand of DNA to search') # use choices array,  use 'plus' and 'minus"
    parser.add_argument('--lcp', type=int, default=12, help='Length of the sequence close to PAM required to be unique')
    parser.add_argument('--dist', type=int, choices=range(1, 6, 1), default=2, help='Minimum Levenshtein edit distance from any other potential target')
    parser.add_argument('--outfile', '-o', type=str, required=True, help='The table of PAM sites and data')
    parser.add_argument('--log', help="Log file", default="guidefinder.log")
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    return parser


def _logger_setup(logfile):
    """Set up logging to a logfile and the terminal standard out.

    Args:
        logfile (str): Log file

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


def main(args=None):
    """Run The complete GuideFinder workflow.

    """
    # Set up logging
    parser = myparser()
    if not args:
        args = parser.parse_args()

    _logger_setup(args.log)
    try:
        pamseq = guidefinder.Pam(args.pamseq, args.pam_orientation)
        if args.tempdir:
            if not os.path.exists(args.tempdir):
                logging.warning("Specified location for tempfile ({}) does not \
                                 exist, using default location.".format(tempdir))
                tempdir = tempfile.mkdtemp(prefix='guidefinder_')

        else:
            tempdir = tempfile.mkdtemp(prefix='guidefinder_', dir=args.tempdir)
            pybedtools.helpers.set_tempdir(tempdir)
        logging.info("Temp directory is: %s", tempdir)
        logging.info("Writing fasta file from genbank file(s)")
        guidefinder.get_fastas(args.gbkfile, tempdir=tempdir)

        logging.info("Identifying PAM sites in the genome")
        possible_targets = []
        input_seqs = SeqIO.parse(os.path.join(tempdir, "forward.fasta"), "fasta")
        for seq in input_seqs:
            possible_targets.append(pamseq.find_targets(seq_obj=seq,
                                                        strand=args.strand,
                                                        target_len=args.targetlength))

        logging.info("Creating a TargetList object")
        targetset = guidefinder.TargetList(targets=possible_targets,
                                             lcp=args.lcp,
                                             levindist=args.dist,
                                             knum=10)
        logging.info("Identifying targets that are unique near the PAM site")
        targetset.find_unique_near_pam()
        logging.info("Index targets for distance search")
        targetset.create_index()
        logging.info("Find targets with the desired minimum edit distance ")
        targetset.get_neighbors()
        logging.info("Convert targets to BED format Pandas Dataframe")
        targets_df = targetset.export_bed()
        logging.info("Parse Genbank features")
        gb_df = guidefinder.get_genbank_features(args.gbkfile)
        logging.info("Identify genome features near PAM sites")
        down, up = guidefinder.get_nearby_feature(targets_df, gb_df)
        # Merge feature columns
        targetfile_columns = list(list(targets_df.values())[0].keys())
        featurefile_columns = list(gb_df.columns) # from genebank feature dataframe
        joined_columns = targetfile_columns + featurefile_columns
        joined_columns.append("distance") # on top of mapping upstream and downstream, we are recording distance in pybed closest. Thus need an extra column
        # merge upstream and downstream output
        logging.info("Merging downstream and upstream features. Columns with \
                      suffix _x represents information from downstream whereas \
                      suffix with _y represent that of upstream")
        merged_down_ups = guidefinder.merge_downstream_upstream(down,
                                                    up,
                                                    joined_columns,
                                                    outputfilename=args.outfile)
        logging.info("Guidefinder completed. Output file is %s" % args.outfile)

    except Exception as e:
        logging.error("Guidefinder terminated with errors. See the log file for details.")
        logging.error(e)
        raise SystemExit(1)
    finally:
        try:
            shutil.rmtree(tempdir)
        except UnboundLocalError:
            pass
        except AttributeError:
            pass
