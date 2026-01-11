import collections
import logging



def _make_feature(seqid, source, type_, start, end, score, strand, phase, attributes):
    attrs = []
    for k, v in attributes.items():
        attrs.append(f"{k}={v}")
    return f"{seqid}\t{source}\t{type_}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{';'.join(attrs)}\n"


def parse_blast_tab6(path, genome_seq, gene_ident=None):
    results = collections.OrderedDict()
    count = 0
    with open(path, 'r') as fh:
        for i, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 12:
                logging.warning(f"Line {i}: unexpected BLAST line with {len(parts)} columns")
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = parts[:12]
            try:
                sstart = int(sstart)
                send = int(send)
            except ValueError:
                logging.warning(f"Line {i}: non-integer coordinates in BLAST sstart/send")
                continue
            start = min(sstart, send)
            end = max(sstart, send)
            strand = '+' if sstart <= send else '-'
            attrs = {
                'ID': f'blast_hit{count}',
                'Target': f'{qseqid} {qstart} {qend}',
                'pident': pident,
                'length': length,
                'evalue': evalue,
                'bitscore': bitscore
            }
            results[f"{start},{end}"] = [strand, '.', 'similarity', attrs]
            count += 1
    return results


def parse_abricate(path, genome_seq, gene_ident=None):
    results = collections.OrderedDict()
    count = 0
    with (open(path, 'r') as fh):
        header = None
        for i, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                header = line.split('\t')
                continue
            if header is None:
                # skip any pre-header content until a header line is encountered
                continue
            parts = line.split('\t')
            if header and len(parts) == len(header):
                row = dict(zip(header, parts))

                try:
                    start = int(row.get('START', '0'))
                    end = int(row.get('END', '0'))
                    strand = row.get('STRAND')
                except ValueError:
                    logging.warning(f"Line {i}: invalid START/END in Abricate line")
                    continue
                seqid = row.get('SEQUENCE')
                gene = row.get('GENE')
                accession =  row.get('ACCESSION') or 'unknown'
                db = row.get('DATABASE')  or 'unknown'
                identity = row.get('%IDENTITY')
                coverage = row.get('%COVERAGE')
                product = row.get('PRODUCT') or 'unkown'
                resistance = row.get('RESISTANCE') or 'unknown'

                attrs = {
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene': gene,
                    'accession': accession,
                    'database': db,
                    'identity': identity,
                    'coverage': coverage,
                    'product': product,
                    'resistance': resistance
                }
                results[f"{start},{end}"] = attrs
                count += 1
            else:
                logging.warning(f"Line {i}: unexpected number of columns in Abricate line")
                continue
    return results


def parse_genemark(path, genome_seq, gene_ident=None):
    results = collections.OrderedDict()
    count = 0
    with open(path, 'r') as fh:
        for i, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                start = int(parts[0])
                stop = int(parts[1])
            except ValueError:
                continue
            strand_tok = parts[2]
            if 'complement' in strand_tok:
                strand = '-'
            else:
                strand = '+'
            attrs = {'ID': f'genemark_hit{count}', 'tool': 'GeneMark'}
            results[f"{start},{stop}"] = [strand, '.', 'CDS', attrs]
            count += 1
    return results


def TabToGFF(input_file, genome_seq, gene_ident='CDS', fmt='blast'):
    # Should be cleaned up to use consistent format names
    fmt = fmt.lower()
    if fmt in ('blast', 'blast_tab6', 'tab6'):
        return parse_blast_tab6(input_file, genome_seq, gene_ident)
    if fmt in ('abricate', 'abricate_tsv', 'abricate_format'):
        return parse_abricate(input_file, genome_seq, gene_ident)
    if fmt in ('genemark', 'gene_mark'):
        return parse_genemark(input_file, genome_seq, gene_ident)
    raise ValueError(f"Unknown format: {fmt}")
