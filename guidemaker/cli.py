"""
GuideMaker: The command line interface
A command line tool to globally design guide RNAs for any CRISPR-Cas system in any small genome
"""

import logging
import argparse
import tempfile
import shutil
import os
import yaml

import pybedtools
from Bio import SeqIO

import guidemaker


def myparser():
    parser = argparse.ArgumentParser(
        description='GuideMaker: globally design guide RNAs for any CRISPR-Cas system in any small genome')
    parser.add_argument('--genbank', '-i', nargs='+', type=str, required=True,
                        help='One or more genbank .gbk  or gzipped .gbk files for a single genome')
    parser.add_argument('--pamseq', '-p', type=str, required=True,
                        help='A short PAM motif to search for, it may use IUPAC ambiguous alphabet')
    parser.add_argument('--outdir', '-o', type=str, required=True,
                        help='The directory for data output')
    parser.add_argument('--pam_orientation', '-r', choices=['5prime', '3prime'], default='5prime',
                        help="PAM position relative to target: 5prime: [PAM][target], 3prime: [target][PAM]. For example, Cas9 is 3prime")
    parser.add_argument('--guidelength', '-l', type=int, default=20, choices=range(10,
                        28, 1), metavar="[10-27]", help='Length of the guide sequence')
    parser.add_argument('--lsr', type=int, default=10, choices=range(0, 28, 1),
                        metavar="[0-27]", help='Length of a seed region near the PAM site required to be unique')
    parser.add_argument('--dist', type=int, choices=range(0, 6, 1),
                        metavar="[0-5]", default=2, help='Minimum hamming distance from any other potential guide')
    parser.add_argument('--before', type=int, default=100, choices=range(1, 501, 1), metavar="[1-500]",
                        help='keep guides this far in front of a feature')
    parser.add_argument('--into', type=int, default=200, choices=range(1, 501, 1), metavar="[1-500]",
                        help='keep guides this far inside (past the start site)of a feature')
    parser.add_argument('--knum', type=int, default=3, choices=range(2, 21, 1),
                        metavar="[2-20]", help='how many sequences similar to the guide to report')
    parser.add_argument('--controls', type=int, default=1000,
                        help='The number or random control RNAs to generate')
    parser.add_argument('--threads', help='The number of cpu threads to use', type=int, default=2)
    parser.add_argument('--log', help="Log file", default="guidemaker.log")
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    parser.add_argument('--restriction_enzyme_list', nargs="*",
                        help='List of sequence representing restriction enzymes', default=[])
    parser.add_argument(
        '--keeptemp', help="Option to keep intermediate files be kept", action='store_true')
    parser.add_argument('--plot', help="Option to genereate guidemaker plots", action='store_true')
    parser.add_argument('--config', help="Path to YAML formatted configuration file, default is " +
                        guidemaker.CONFIG_PATH, default=guidemaker.CONFIG_PATH)
    parser.add_argument('-V', '--version', action='version',
                        version="%(prog)s (" + guidemaker.__version__ + ")")
    return parser


def parserval(args):
    assert(args.lsr <= args.guidelength), "The length of sequence near the PAM .i.e seed sequence that must be less than the guide length"
    # Campylobacter jejuni Cas9 (CjCas9) has a 8bp long 5’-NNNNRYAC-3’ PAM site
    assert(1 < len(args.pamseq) < 9), "The length of the PAM sequence must be between 2-8"


def _logger_setup(logfile):
    """
    Set up logging to a logfile and the terminal standard out.

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


def main(arglist: list = None):
    """Run The complete GuideMaker workflow."""
    # Set up logging
    parser = myparser()
    args = parser.parse_args(arglist)
    parserval(args)

    _logger_setup(args.log)

    try:
        with open(args.config) as cf:
            config = yaml.safe_load(cf)
    except:
        print("Could not parse the configuration file.")
        raise SystemExit(1)

    try:
        logging.info("Configuration data loaded from {}:".format(args.config))
        logging.info(config)
    except:
        print("Could find config file, exiting.")
        raise SystemExit(1)

    try:

        if args.tempdir:
            if not os.path.exists(args.tempdir):
                logging.warning("Specified location for tempfile (%s) does not \
                                 exist, using default location." % args.tempdir)
                tempdir = tempfile.mkdtemp(prefix='guidemaker_')
        else:
            tempdir = tempfile.mkdtemp(prefix='guidemaker_', dir=args.tempdir)
            pybedtools.helpers.set_tempdir(tempdir)
        logging.info("Temp directory is: %s" % (tempdir))
        logging.info("Writing fasta file from genbank file(s)")
        fastapath = guidemaker.get_fastas(args.genbank, tempdir=tempdir)
        logging.info("Identifying PAM sites in the genome")
        pamobj = guidemaker.core.PamTarget(args.pamseq, args.pam_orientation)
        seq_record_iter = SeqIO.parse(fastapath, "fasta")
        pamtargets = pamobj.find_targets(
            seq_record_iter=seq_record_iter, target_len=args.guidelength)
        tl = guidemaker.core.TargetProcessor(
            targets=pamtargets, lsr=args.lsr, hammingdist=args.dist, knum=args.knum)
        lengthoftl = len(tl.targets)
        logging.info("Checking guides for restriction enzymes")
        tl.check_restriction_enzymes(restriction_enzyme_list=args.restriction_enzyme_list)
        logging.info("Number of guides removed after checking for restriction enzymes: %d",
                     (lengthoftl - len(tl.targets)))
        logging.info("Identifing guides that are unique near the PAM site")
        tl.find_unique_near_pam()
        logging.info("Number of guides with non unique seed sequence: %d",
                     (tl.targets.isseedduplicated.sum()))
        logging.info("Indexing all potential guide sites: %s. This is the longest step." %
                     len(list(set(tl.targets['target'].tolist()))))
        tl.create_index(num_threads=args.threads, configpath=args.config)
        logging.info(
            "Identifying guides that have a hamming distance <= %s to all other potential guides", str(args.dist))
        tl.get_neighbors(num_threads=args.threads, configpath=args.config)
        logging.info("Formatting data for BedTools")
        tf_df = tl.export_bed()
        logging.info("Create GuideMaker Annotation object")
        anno = guidemaker.core.Annotation(genbank_list=args.genbank,
                                           target_bed_df=tf_df)
        logging.info("Identify genomic features")
        anno._get_genbank_features()
        logging.info("Total number of CDS/locus in the input genome: %d" % anno.locuslen())
        logging.info("Find genomic features closest the guides")
        anno._get_nearby_features()
        logging.info("Select guides that start between +%s and -%s of a feature start" %
                     (args.before, args.into))
        anno._filter_features(before_feat=args.before, after_feat=args.into)
        logging.info("Select description columns")
        anno._get_qualifiers(configpath=args.config)
        logging.info("Format the output")
        prettydf = anno._format_guide_table(tl)
        fd_zero = prettydf['Feature distance'].isin([0]).sum()
        logging.info("Number of Guides within a gene coordinates i.e. zero Feature distance: %d", fd_zero)
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        csvpath = os.path.join(args.outdir, "targets.csv")
        prettydf.to_csv(csvpath, index=False)
        logging.info("creating random control guides")
        contpath = os.path.join(args.outdir, "controls.csv")
        seq_record_iter = SeqIO.parse(fastapath, "fasta")
        cmin, cmed, randomdf = tl.get_control_seqs(
            seq_record_iter, configpath=args.config, length=args.guidelength, n=args.controls, num_threads=args.threads)
        logging.info("Number of random control searched: %d" % tl.ncontrolsearched)
        logging.info("Percentage of GC content in the input genome: " +
                     "{:.2f}".format(tl.gc_percent))
        logging.info("Total length of the genome: " + "{:.1f} MB".format(tl.genomesize))
        logging.info("created %i control guides with a minimum distance of %d and a median distance of %d" % (
            args.controls, cmin, cmed))
        randomdf.to_csv(contpath)
        logging.info("guidemaker completed, results are at %s" % args.outdir)
        logging.info("PAM sequence: %s" % args.pamseq)
        logging.info("PAM orientation: %s" % args.pam_orientation)
        logging.info("Genome strand(s) searched: %s" % "both")
        logging.info("Total PAM sites considered: %d" % lengthoftl)
        logging.info("Guide RNA candidates found: %d" % len(prettydf))
    except Exception as e:
        logging.error("GuideMaker terminated with errors. See the log file for details.")
        logging.error(e)
        raise SystemExit(1)
    try:
        if args.plot:
            logging.info("Creating Plots...")
            guidemaker.core.GuideMakerPlot(prettydf=prettydf, outdir=args.outdir)
            logging.info("Plots saved at: %s" % args.outdir)
    except Exception as e:
        logging.error(e)
        raise SystemExit(1)
    try:
        if not args.keeptemp:
            shutil.rmtree(tempdir)
    except UnboundLocalError as e:
        logging.error(e)
        raise SystemExit(1)
    except AttributeError as e:
        logging.error(e)
        raise SystemExit(1)
