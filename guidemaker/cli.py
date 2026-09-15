"""
GuideMaker: The command line interface
A command line Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems
"""

import logging

import argparse
import tempfile
import shutil
import os
import yaml
import textwrap


import pybedtools
from Bio import SeqIO

import guidemaker


def myparser():
    parser = argparse.ArgumentParser(
        description='GuideMaker: Software to design gRNAs pools in non-model genomes and CRISPR-Cas systems',
        epilog=textwrap.dedent(''' To run the web app locally, in terminal run:
        -----------------------------------------------------------------------
        streamlit run ''' + str(guidemaker.WEB_APP) + '''
        -----------------------------------------------------------------------'''))
    parser.add_argument('--genbank', '-i', nargs='+', type=str, required=False,
                        help='One or more genbank .gbk  or gzipped .gbk files for a single genome. Provide this or GFF/GTF and fasta files')
    parser.add_argument('--fasta', '-f', nargs='+', type=str, required=False,
                        help='One or more fasta or gzipped fasta files for a single genome. If using a fasta, a GFF/GTF file must also be provided but not a genbank file.')
    parser.add_argument('--gff', '-g', nargs='+', type=str, required=False,
                        help='One or more GFF or GTF files (optionally gzipped) for a single genome. If using a GFF/GTF a fasta file  must also be provided but not a genbank file.')
    parser.add_argument('--pamseq', '-p', type=str, required=True,
                        help='A short PAM motif to search for, it may use IUPAC ambiguous alphabet')
    parser.add_argument('--outdir', '-o', type=str, required=True,
                        help='The directory for data output')
    parser.add_argument('--raw_output_only', action='store_true',
                        help='if selected only the raw guide RNAs their positions will be returned the meet lsr and dist criteria')
    parser.add_argument('--pam_orientation', '-r', choices=['5prime', '3prime'], default='3prime',
                        help="The PAM position relative to the target: 5prime: [PAM][target], 3prime: [target][PAM]. For example, SpCas9 is 3prime. Default: '3prime'.")
    parser.add_argument('--guidelength', '-l', type=int, default=20, choices=range(10,
                        28, 1), metavar="[10-27]", help='Length of the guide sequence. Default: 20.')
    parser.add_argument('--lsr', type=int, default=10, choices=range(0, 28, 1),
                        metavar="[0-27]", help='Length of a seed region near the PAM site required to be unique. Default: 10.')
    parser.add_argument('--dtype', type=str, choices=['hamming', 'leven'], default='hamming',
                        help='Select the distance type. Default: hamming.')
    parser.add_argument('--dist', type=int, choices=range(0, 6, 1),
                        metavar="[0-5]", default=2, help='Minimum edit distance from any other potential guide. Default: 2.')
    parser.add_argument('--before', type=int, default=100, choices=range(1, 501, 1), metavar="[1-500]",
                        help='keep guides this far in front of a feature. Default: 100.')
    parser.add_argument('--into', type=int, default=200, metavar="[1+]",
                        help='keep guides this far inside (past the start site) of a feature. Default: 200.')
    parser.add_argument('--knum', type=int, default=5, choices=range(2, 21, 1),
                        metavar="[2-20]", help='how many sequences similar to the guide to report. Default: 5.')
    parser.add_argument('--controls', type=int, default=1000, choices=range(0, 100001, 1),  metavar="[0-100000]",
                        help='Number of random control RNAs to generate. Default: 1000.')
    parser.add_argument('--threads', help='The number of cpu threads to use. Default: 2', type=int, default=2)
    parser.add_argument('--log', help="Log file", default="guidemaker.log")
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    parser.add_argument('--restriction_enzyme_list', nargs="*",
                        help='List of sequence representing restriction enzymes. Default: None.', default=[])
    parser.add_argument('--attribute_key', type=str,
                        help='the attribute key in column 9 of the GFF/GTF file to use for filtering. Default: ID', default="ID")
    parser.add_argument('--filter_by_attribute', nargs="*",
                        help='List of locus ids. Default: None.', default=[])
    parser.add_argument('--doench_efficiency_score',help="On-target scoring from Doench et al. 2016 - only for NGG PAM and guidelength=25: Default: None.", action='store_true')
    parser.add_argument('--cfd_score',help='CFD score for assessing off-target activity of gRNAs with NGG pam: Default: None.', action='store_true')
    parser.add_argument('--keeptemp', help="Option to keep intermediate files be kept", action='store_true')
    parser.add_argument('--plot', help="Option to create GuideMaker plots", action='store_true')
    parser.add_argument('--config', help="Path to YAML formatted configuration file, default is " +
                        str(guidemaker.CONFIG_PATH), default=str(guidemaker.CONFIG_PATH))
    parser.add_argument('-V', '--version', action='version',
                        version="%(prog)s (" + guidemaker.__version__ + ")")
    return parser



def parserval(args):
    try:
        assert(args.lsr <= args.guidelength), "The length of sequence near the PAM .i.e seed sequence that must be less than the guide length"
        # Campylobacter jejuni Cas9 (CjCas9) has a 8bp long 5’-NNNNRYAC-3’ PAM site
        assert(1 < len(args.pamseq) < 9), "The length of the PAM sequence must be between 2-8"
        assert ((args.genbank is not None and args.fasta is None and args.gff is None) or
                (args.genbank is None and args.fasta is not  None and args.gff is not None) or
                ((args.genbank is not None or args.fasta is not None) and args.raw_output_only)),"Please provide either Genbank files or Fasta and GFF files. If raw_output_onl is selected Genbank or Fasta files are requried."
    except AssertionError as err:
        raise err

def _logger_setup(logfile):
    """
    Set up logging to a logfile and the terminal standard out.

    Args:
        logfile (str): Log file
    """
    try:
        # create logger
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        # create console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create file handler

        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
        # tell the handler to use this format
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(fh)
        logger.addHandler(ch)
        return logger
    except Exception as e:
        print("An error occurred setting up logging")
        raise e


def main(arglist: list = None):
    """Run The complete GuideMaker workflow."""
    # Set up logging
    parser = myparser()
    args = parser.parse_args(arglist)
    logger = _logger_setup(args.log)
    parserval(args)



    try:
        with open(args.config) as cf:
            config = yaml.safe_load(cf)
    except:
        print("Could not parse the configuration file.")
        raise SystemExit(1)

    try:
        logger.info("Configuration data loaded from {}:".format(args.config))
        logger.info(config)
    except:
        print("Could not find config file, exiting.")
        raise SystemExit(1)

    try:

        if args.tempdir:
            if not os.path.exists(args.tempdir):
                logger.warning("Specified location for tempfile (%s) does not \
                                 exist, using default location." % args.tempdir)
                os.mkdir(args.tempdir)
                tempdir = args.tempdir
            else:
                tempdir = tempfile.mkdtemp()
        else:
            tempdir = tempfile.mkdtemp(prefix='guidemaker_', dir=args.tempdir)
            pybedtools.helpers.set_tempdir(tempdir)
        logger.info("Temp directory is: %s" % (tempdir))
        if args.genbank:
            logger.info("Writing fasta file from genbank file(s)")
            fastapath = guidemaker.get_fastas(args.genbank, input_format="genbank", tempdir=tempdir)
        elif args.fasta:
            fastapath = guidemaker.get_fastas(args.fasta, input_format="fasta", tempdir=tempdir)
        logger.info("Identifying PAM sites in the genome")
        pamobj = guidemaker.core.PamTarget(args.pamseq, args.pam_orientation, args.dtype)
        seq_record_iter = SeqIO.parse(fastapath, "fasta")
        pamtargets = pamobj.find_targets(
            seq_record_iter=seq_record_iter, target_len=args.guidelength)
        tl = guidemaker.core.TargetProcessor(
            targets=pamtargets, lsr=args.lsr, editdist=args.dist, knum=args.knum)
        lengthoftl = len(tl.targets)
        logger.info("Checking guides for restriction enzymes")
        tl.check_restriction_enzymes(restriction_enzyme_list=args.restriction_enzyme_list)
        logger.info("Number of guides removed after checking for restriction enzymes: %d",
                     (lengthoftl - len(tl.targets)))
        logger.info("Identifing guides that are unique near the PAM site")
        tl.find_unique_near_pam()
        logger.info("Number of guides with non unique seed sequence: %d",
                     (tl.targets.isseedduplicated.sum()))
        logger.info("Indexing all potential guide sites: %s. This is the longest step." %
                     len(list(set(tl.targets['target'].tolist()))))
        tl.create_index(num_threads=args.threads, configpath=args.config)
        logger.info(
            "Identifying guides that have a hamming distance <= %s to all other potential guides", str(args.dist))
        tl.get_neighbors(num_threads=args.threads, configpath=args.config)
        logger.info("Formatting data for BedTools")
        tf_df = tl.export_bed()
        if args.raw_output_only:
            tf_df.to_csv(os.path.join(args.outdir, "rawguides.csv.gz"), index=False, header=["Chromosome", "Start", "Stop","gRNA", "Strand"])
            logger.info("Raw guides options was selected Guidemaker, so has completed opperations")
            raise SystemExit(0)

        logger.info("Create GuideMaker Annotation object")
        if args.genbank:
            anno = guidemaker.core.Annotation(annotation_list=args.genbank, annotation_type="genbank",
                                              target_bed_df=tf_df)
        elif args.gff:
            anno = guidemaker.core.Annotation(annotation_list=args.gff, annotation_type="gff",
                                              target_bed_df=tf_df)
        logger.info("Identify genomic features")
        anno.get_annotation_features()
        logger.info("Total number of %s in the input genome: %d" % anno.locuslen())
        logger.info("Find genomic features closest the guides")
        anno._get_nearby_features()
        logger.info("Select guides that start between +%s and -%s of a feature start" %
                     (args.before, args.into))
        anno._filter_features(before_feat=args.before, after_feat=args.into)
        logger.info("Select description columns")
        anno._get_qualifiers(configpath=args.config)
        logger.info("Format the output")
        anno._format_guide_table(tl)
        prettydf = anno._filterlocus(args.attribute_key, args.filter_by_attribute)
        # prettydf = anno.filter_pretty_df
        if args.doench_efficiency_score:
            logger.info("Creating Efficiency Score based on Doench et al. 2016 - only for NGG PAM...")
            prettydf = guidemaker.core.get_doench_efficiency_score(df=prettydf, pam_orientation=args.pam_orientation, num_threads=args.threads)

        if args.cfd_score:
            logger.info("Calculating CFD score for assessing off-target activity of gRNAs")
            prettydf = guidemaker.core.cfd_score(df=prettydf)

        fd_zero = prettydf['Feature distance'].isin([0]).sum()
        logger.info("Number of Guides within a gene coordinates i.e. zero Feature distance: %d", fd_zero)
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        csvpath = os.path.join(args.outdir, "targets.csv.gz")
        prettydf.to_csv(csvpath, index=False)
        if args.controls > 0:
            logger.info("Creating random control guides")
            contpath = os.path.join(args.outdir, "controls.csv.gz")
            seq_record_iter = SeqIO.parse(fastapath, "fasta")
            cmin, cmed, randomdf = tl.get_control_seqs(seq_record_iter,
                                                       configpath=args.config,
                                                       length=args.guidelength,
                                                       n=args.controls,
                                                       num_threads=args.threads)
            randomdf.to_csv(contpath)
            logger.info("Number of random control searched: %d" % tl.ncontrolsearched)
            logger.info("Created %i control guides with a minimum distance of %d and a median distance of %d" % (
                args.controls, cmin, cmed))
            logger.info("Percentage of GC content in the input genome: " +
                        "{:.2f}".format(tl.gc_percent))
            logger.info("Total length of the genome: " + "{:.1f} MB".format(tl.genomesize))


        logger.info("GuideMaker completed, results are at %s" % args.outdir)
        logger.info("PAM sequence: %s" % args.pamseq)
        logger.info("PAM orientation: %s" % args.pam_orientation)
        logger.info("Genome strand(s) searched: %s" % "both")
        logger.info("Total PAM sites considered: %d" % lengthoftl)
        logger.info("Guide RNA candidates found: %d" % len(prettydf))
    except Exception as e:
        logger.exception("GuideMaker terminated with errors. See the log file for details.")
        raise SystemExit(1)
    try:
        if args.plot:
            logger.info("Creating Plots...")
            guidemaker.core.GuideMakerPlot(prettydf=prettydf, outdir=args.outdir)
            logger.info("Plots saved at: %s" % args.outdir)
    except Exception as e:
        logger.exception(e)
        raise SystemExit(1)
    try:
        if not args.keeptemp:
            shutil.rmtree(tempdir)
    except UnboundLocalError as e:
        logger.exception(e)
        raise SystemExit(1)
    except AttributeError as e:
        logger.exception(e)
        raise SystemExit(1)
