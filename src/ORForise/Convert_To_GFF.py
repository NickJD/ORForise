import argparse
import logging
from datetime import datetime
import os
import sys

try:
    from utils import *
    from ORForise.src.ORForise.Aux.TabToGFF import TabToGFF
except ImportError:
    from ORForise.utils import *
    from ORForise.Aux.TabToGFF.TabToGFF import TabToGFF


def setup_logging(outdir, verbose=False):
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    logfile = None
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    # clear existing handlers to avoid duplicates when running repeatedly
    logger.handlers = []
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # Only create a file handler (and thus the logfile) when verbose is enabled
    if verbose:
        logfile = os.path.join(outdir, f'convert_to_gff_{ts}.log')
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    # Always add a stdout handler
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG if verbose else logging.INFO)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    return logfile


def write_gff(outpath, genome_ID, genome_DNA, input_annotation, fmt, features):
    with open(outpath, 'w') as out:
        out.write('##gff-version\t3\n')
        out.write('#\tConvert_To_GFF\n')
        out.write('#\tRun Date: ' + str(datetime.now()) + '\n')
        # Only include genome DNA line if a path was provided
        if genome_DNA:
            out.write('##Genome DNA File:' + genome_DNA + '\n')
        out.write('##Original File: ' + input_annotation + '\n')
        for pos, data in features.items():
            pos_ = pos.split(',')
            start = pos_[0]
            stop = pos_[-1]
            strand = data['strand']
            if fmt == 'abricate': # Currently only supports abricate format
                info = 'abricate_anotation;accession='+data['accession']+';database='+data['database']+';identity='+str(data['identity'])+';coverage='+str(data['coverage'])+';product='+data['product']+';resistance='+data['resistance']
            entry = f"{data['seqid']}\t{fmt}\t{'CDS'}\t{start}\t{stop}\t.\t{strand}\t.\t{'ID='}{info}\n"
            out.write(entry)


def load_genome(genome_fasta):
    genome_seq = ''
    genome_ID = 'unknown'
    with open(genome_fasta, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                genome_ID = line.split()[0].lstrip('>')
            else:
                genome_seq += line
    return genome_ID, genome_seq


def main():
    print(WELCOME)

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': Convert-To-GFF Run Parameters')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')

    required.add_argument('-i', dest='input_annotation', required=True, help='Input annotation file (tabular)')
    required.add_argument('-fmt', dest='format', required=True, help='Input format: blast, abricate, genemark')
    required.add_argument('-o', dest='output_dir', required=True, help='Output directory')

    optional = parser.add_argument_group('Optional Arguments')
    # Make genome DNA optional: if not provided we operate without genome sequence
    required.add_argument('-dna', dest='genome_DNA', required=False, help='Genome DNA file (.fa)')
    optional.add_argument('-gi', dest='gene_ident', default='CDS', required=False, help='Gene identifier types to extract (unused)')
    optional.add_argument('--verbose', dest='verbose', action='store_true', help='Verbose logging with logfile')

    options = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    logfile = setup_logging(options.output_dir, verbose=options.verbose)
    logging.info('Starting Convert_To_GFF')
    # Log genome DNA only if provided
    if options.genome_DNA:
        logging.info('Genome DNA: %s', options.genome_DNA)
    else:
        logging.info('Genome DNA: (not provided)')
    logging.info('Input annotation: %s', options.input_annotation)
    logging.info('Format: %s', options.format)

    # If a genome fasta was provided, load it; otherwise proceed without genome sequence
    if options.genome_DNA:
        if not os.path.exists(options.genome_DNA):
            logging.error('Genome DNA file does not exist: %s', options.genome_DNA)
            sys.exit(1)
        genome_ID, genome_seq = load_genome(options.genome_DNA)
    else:
        # Derive a sensible genome_ID from the annotation filename and leave sequence empty
        genome_ID = os.path.splitext(os.path.basename(options.input_annotation))[0]
        genome_seq = ''

    try:
        # Build genome map expected by TabToGFF: mapping genome_ID -> tuple(sequence, ...)
        genome_map = {genome_ID: (genome_seq,)}
        features = TabToGFF(options.input_annotation, genome_map, options.gene_ident, fmt=options.format)
    except Exception as e:
        logging.exception('Error parsing input annotation')
        sys.exit(1)

    #features = sortORFs(features) - Not sorting for now to preserve original order
    basename = os.path.basename(options.input_annotation)
    dot = basename.rfind('.')
    if dot != -1:
        outname = basename[:dot] + '.gff'
    else:
        outname = basename + '.gff'
    outgff = os.path.join(options.output_dir, outname)
    # Pass the original genome path if provided, else pass None so headers adapt
    genome_DNA_path = options.genome_DNA if options.genome_DNA else None
    write_gff(outgff, genome_ID, genome_DNA_path, options.input_annotation, options.format, features)
    logging.info('Wrote GFF to %s', outgff)
    logging.info('Logfile: %s', logfile)

if __name__ == '__main__':
    main()
