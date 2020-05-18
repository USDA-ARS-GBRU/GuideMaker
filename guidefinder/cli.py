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
    parser = argparse.ArgumentParser(description='GuideFinder: globally design guide RNAs for any CRISPR-Cas system in any small genome')
    parser.add_argument('--genbank', '-i', nargs='+', type=str, required=True, help='One or more genbank .gbk  or gzipped .gbk files for a single genome')
    parser.add_argument('--pamseq', '-p', type=str, required=True, help='A short PAM motif to search for, it may use IUPAC ambiguous alphabet')
    parser.add_argument('--outfile', '-o', type=str, required=True, help='The table of PAM sites and data')
    parser.add_argument('--pam_orientation', '-r', choices=['5prime', '3prime'], default='5prime', help="PAM position relative to target: 5prime: [PAM][target], 3prime: [target][PAM], Cas9 for example is 5prime")
    parser.add_argument('--guidelength', '-l', type=int, default=20, help='Length of the guide sequence')
    parser.add_argument('--strand', '-s', choices=['forward','reverse', 'both'], default='both', help='Strand of DNA to search')
    parser.add_argument('--lcp', type=int, default=10, help='Length of the guide closest to  the PAM required to be unique')
    parser.add_argument('--dist', type=int, choices=range(1, 6, 1), default=2, help='Minimum hamming distance from any other potential guide')
    parser.add_argument('--before', type=int, default=100,
                        help='keep guides this far in front of a feature')
    parser.add_argument('--into', type=int, default=200,
                        help='keep guides this far inside (past the start site)of a feature')
    parser.add_argument('--threads', help='The number of cpu threads to use', type=int, default=2)
    parser.add_argument('--log', help="Log file", default="guidefinder.log")
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s (" + guidefinder.__version__ + ")")
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

        if args.tempdir:
            if not os.path.exists(args.tempdir):
                logging.warning("Specified location for tempfile (%s) does not \
                                 exist, using default location." % tempdir )
                tempdir = tempfile.mkdtemp(prefix='guidefinder_')
        else:
            tempdir = tempfile.mkdtemp(prefix='guidefinder_', dir=args.tempdir)
            pybedtools.helpers.set_tempdir(tempdir)
        logging.info("Temp directory is: %s" % (tempdir))
        logging.info("Writing fasta file from genbank file(s)")

        fastapath = guidefinder.get_fastas(args.genbank, tempdir=tempdir)
        logging.info("Identifying PAM sites in the genome")
        pamobj = guidefinder.core.Pam(args.pamseq, args.pam_orientation)
        fastahandle = SeqIO.read(fastapath, "fasta")
        pamtargets = pamobj.find_targets(seqrecord_obj=fastahandle, strand="both", target_len=args.guidelength)
        tl = guidefinder.core.TargetList(targets=pamtargets, lcp=args.lcp, hammingdist=args.dist, knum=2)
        logging.info("Identifing guides that are unique near the PAM site")
        tl.find_unique_near_pam()
        logging.info("Indexing all potential guide sites: %s. This is the longest step." % len(tl.targets))
        tl.create_index(num_threads=args.threads)
        logging.info("Identifying guides that have a hamming distance <= %s to all other potential guides", str(args.dist))
        tl.get_neighbors(num_threads=args.threads)
        logging.info("Formatting data for BedTools")
        tf_df = tl.export_bed()
        logging.info("Create Guidefinder Annotation object")
        anno = guidefinder.core.Annotation(genbank_list=args.genbank,
                                           target_bed_df=tf_df)
        logging.info("Identify genomic features")
        anno._get_genbank_features()
        logging.info("Find genomic features closest the guides")
        anno._get_nearby_features()
        logging.info("Select guides that start between +%s and -%s of a feature start" % (args.before, args.into))
        anno._filter_features(before_feat=args.before, after_feat=args.into)
        logging.info("Select description columns")
        anno._get_qualifiers()
        logging.info("Format the output")
        prettydf = anno._format_guide_table(tl)
        prettydf.to_csv(args.outfile)
        logging.info("Guidefinder completed, results are at %s" % args.outfile)
        logging.info("PAM sequence: %s" % args.pamseq)
        logging.info("PAM orientation: %s" % args.pam_orientation)
        logging.info("Genome strand(s) searched: %s" % args.strand)
        logging.info("Total PAM sites considered: %d" % len(tl.targets))
        logging.info("Guide RNA candidates found: %d" % len(prettydf))

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
