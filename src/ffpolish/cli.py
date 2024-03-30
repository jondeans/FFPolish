"""FFPolish CLI Entry Point."""

from pathlib import Path
import importlib.metadata

from argparse import ArgumentParser

from ffpolish.filter import filter
from ffpolish.extract import extract

from enum import Enum

__version__ = importlib.metadata.version(__package__ or __name__)


def main():

    parser = ArgumentParser(description="FFPolish - Filter Artifacts From FFPE Variant Calls")
    subparsers = parser.add_subparsers(dest="command")

    parser.add_argument("--loglevel", type=str, default="info", choices=["debug", "info", "warning", "error", "critical"], help="Set the logging level")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-o", "--outdir", default="results", help='Output directory')
    parser.add_argument("-p", "--prefix", default=None, help='Output prefix (default: basename of BAM)')

    filter_parser = subparsers.add_parser('filter', help='Filter variants')
    filter_parser.add_argument('--retrain', action="store_true", help='Retrain model with new data (in hdf5 format)')
    filter_parser.add_argument('--grid-search', action='store_true', default=False, help='Perform grid search when retraining model')
    filter_parser.add_argument('-c', '--cores', default=1, help='Number of cores to use for grid search (default: %(default)s)')
    filter_parser.add_argument('-s', '--seed', default=24159, help="Seed for retraining (default: %(default)s)")
    filter_parser.add_argument('REF', help="Reference genome FASTA.")
    filter_parser.add_argument('VCF', help="VCF to filter.")
    filter_parser.add_argument('BAM', help="Tumor BAM file.")

    extract_parser = subparsers.add_parser('extract', help='Extract features for re-training')
    extract_parser.add_argument('--skip-bam-readcount', action='store_true', default=False, help='Skip bam_readcount on sample')
    extract_parser.add_argument('--labels', default=None, help='BED file of true variants (to create pickle file of true variants)')
    extract_parser.add_argument('--pkl', default=None, help='Pickle file of true variants')
    extract_parser.add_argument('REF', help="Reference genome FASTA.")
    extract_parser.add_argument('VCF', help="VCF to filter.")
    extract_parser.add_argument('BAM', help="Tumor BAM file.")

    args = parser.parse_args()

    if "REF" in args and isinstance(args.REF, str):
        args.REF = Path(args.REF)
    if "VCF" in args and isinstance(args.VCF, str):
        args.VCF = Path(args.VCF)
    if "BAM" in args and isinstance(args.BAM, str):
        args.BAM = Path(args.BAM)
    if "outdir" in args and isinstance(args.outdir, str):
        args.outdir = Path(args.outdir)
    if "prefix" in args and isinstance(args.prefix, str):
        args.prefix = Path(args.prefix)

    if "retrain" in args:
        args.retrain = bool(args.retrain)
    if "grid_search" in args:
        args.grid_search = bool(args.grid_search)
    if "cores" in args:
        args.cores = int(args.cores)
    if "seed" in args:
        args.seed = int(args.seed)

    if args.command == 'filter':
        filter(args.REF, args.VCF, args.BAM, args.outdir, args.prefix, args.retrain, args.grid_search, args.cores, args.seed)
    elif args.command == 'extract':
        extract(args.ref, args.bam, args.outdir, args.prefix, args.labels)


if __name__ == '__main__':
    main()
