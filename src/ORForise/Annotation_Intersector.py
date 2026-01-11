from importlib import import_module
import argparse
from collections import OrderedDict
from datetime import date
import sys, gzip
import os
import logging

# Ensure logging prints to stdout by default so info/debug messages are visible when running the script
if not logging.getLogger().handlers:
    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(levelname)s: %(message)s')

try:
    from utils import *
except ImportError:
    from .utils import *

################################


def gff_writer(genome_ID, genome_DNA, reference_annotation, reference_tool, ref_gene_set, additional_annotation, additional_tool, genes_To_Keep_by_contig, output_file, gene_ident=None):
    # genes_To_Keep_by_contig: {contig: {pos: data}}
    # Expand user (~) and ensure output directory exists
    output_file = os.path.expanduser(output_file)
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Always open the file and write the header first. Use a broad try so we can log any issue.
    try:
        with open(output_file, 'w', encoding='utf-8') as write_out:
            write_out.write("##gff-version\t3\n#\tAnnotation-Intersector\n#\tRun Date:" + str(date.today()) + '\n')
            if genome_DNA:
                write_out.write("##Genome DNA File:" + genome_DNA + '\n')
            write_out.write("##Original File: " + reference_annotation + "\n##Intersecting File: " + additional_annotation + '\n')

            entries_written = 0

            # If genes_To_Keep_by_contig is falsy or empty, attempt to write reference features directly as fallback
            if not genes_To_Keep_by_contig or all(len(v) == 0 for v in genes_To_Keep_by_contig.values()):
                write_out.write(f"# No kept genes to write (0 entries).\n")
                write_out.write(f"# Falling back to writing reference features with coverage=0.\n")

                # Parse reference annotation and write features matching gene_ident
                try:
                    if reference_annotation.endswith('.gz'):
                        rf = gzip.open(reference_annotation, 'rt')
                    else:
                        rf = open(reference_annotation, 'r', encoding='unicode_escape')
                    with rf:
                        for line in rf:
                            line = line.rstrip('\n')
                            if not line or line.startswith('#'):
                                continue
                            parts = line.split('\t')
                            if len(parts) < 9:
                                continue
                            seqid = parts[0]
                            ftype = parts[2]
                            try:
                                gene_types = gene_ident.split(',') if gene_ident else ['CDS']
                            except Exception:
                                gene_types = ['CDS']
                            if ftype not in gene_types and not ('CDS' in gene_types and ftype == 'CDS'):
                                continue
                            try:
                                start = parts[3]
                                stop = parts[4]
                                strand = parts[6]
                                info = parts[8]
                            except Exception:
                                continue
                            # write entry with coverage 0 and empty additional annotation
                            entry = f"{seqid}\t{os.path.splitext(os.path.basename(reference_annotation))[0]}\t{ftype}\t{start}\t{stop}\t.\t{strand}\t.\tID=Original_Annotation={info};Additional_Annotation=;Coverage=0\n"
                            write_out.write(entry)
                            entries_written += 1
                except Exception as e:
                    logging.warning('Fallback parse of reference annotation failed: %s', e)

                write_out.flush()
                logging.info('Wrote %d fallback reference entries to %s', entries_written, output_file)
                return

            for contig, genes in genes_To_Keep_by_contig.items():
                # Use basename without extension for the source field
                ref = os.path.splitext(os.path.basename(reference_annotation))[0].split('_')[0]
                for pos, data in genes.items():
                    try:
                        pos_ = pos.split(',')
                        start = pos_[0]
                        stop = pos_[-1]
                        strand = data[0]
                        # Ensure indices exist and are strings
                        add_ann = str(data[4]) if len(data) > 4 else ''
                        orig_ann = str(data[5]) if len(data) > 5 else ''
                        entry = (
                            contig + '\t' + ref + '\t' + data[2] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation=' + orig_ann + ';Additional_Annotation=' + add_ann + ';Coverage=' + str(
                                data[1]) + '\n')
                        write_out.write(entry)
                        entries_written += 1
                    except Exception as e:
                        # Log the bad entry and continue
                        logging.warning('Skipping bad GFF entry for contig %s pos %s: %s', contig, pos, e)
                        continue

            write_out.flush()
            logging.info('Wrote %d GFF entries to %s', entries_written, output_file)
    except OSError as e:
        logging.error("Cannot write to output file %s: %s", output_file, e)
        sys.exit(1)


def _get_opt(options, *names):

    for n in names:
        if hasattr(options, n):
            return getattr(options, n)
    return None


def _parse_pos(pos_str):
    try:
        s, e = pos_str.split(',')
        return int(s), int(e)
    except Exception:
        return None, None


def _write_discordance_report(report_path, entries):
    # Summarise discordance entries instead of writing each row (GFFs keep the full detail)
    report_path = os.path.expanduser(report_path)
    out_dir = os.path.dirname(report_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    try:
        from collections import Counter
        total = len(entries) if entries is not None else 0
        status_counts = Counter()
        contig_counts = Counter()
        cov_values = []
        for e in (entries or []):
            st = str(e.get('status', 'unknown'))
            status_counts[st] += 1
            contig = str(e.get('contig', ''))
            if contig:
                contig_counts[contig] += 1
            # attempt to extract numeric coverage if present
            try:
                cov = float(e.get('coverage', 0) if e.get('coverage', '') != '' else 0)
                cov_values.append(cov)
            except Exception:
                try:
                    # sometimes coverage might be a string like '12.34' or '0.00'
                    cov_values.append(float(str(e.get('coverage', '0')).strip()))
                except Exception:
                    pass

        avg_cov = (sum(cov_values) / len(cov_values)) if cov_values else 0.0
        nonzero_covs = [v for v in cov_values if v != 0]
        nonzero_avg_cov = (sum(nonzero_covs) / len(nonzero_covs)) if nonzero_covs else 0.0

        with open(report_path, 'w', encoding='utf-8') as fh:
            fh.write('metric\tvalue\n')
            fh.write(f'total_entries\t{total}\n')
            fh.write(f'unique_contigs_with_discordance\t{len(contig_counts)}\n')
            fh.write(f'average_coverage_reported\t{avg_cov:.2f}\n')
            fh.write(f'non-zero_average_coverage_reported\t{nonzero_avg_cov:.2f}\n')
            # statuses
            for status, cnt in status_counts.most_common():
                fh.write(f'status::{status}\t{cnt}\n')
            # per-contig counts (only top 50 to avoid huge files)
            fh.write('# per-contig discordance counts (top 50)\n')
            for contig, cnt in contig_counts.most_common(50):
                fh.write(f'contig::{contig}\t{cnt}\n')
        logging.info('Wrote discordance summary report: %s', report_path)
    except OSError as e:
        logging.error('Cannot write discordance report %s: %s', report_path, e)
        sys.exit(1)


def _write_discordance_gff(report_path, entries, reference_annotation_basename):
    """Write a list of discordance entries (dicts) to a GFF file."""
    report_path = os.path.expanduser(report_path)
    out_dir = os.path.dirname(report_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    try:
        with open(report_path, 'w', encoding='utf-8') as fh:
            fh.write('##gff-version\t3\n')
            fh.write('#\tAnnotation-Intersector discordance report\n')
            fh.write('#\tRun Date:' + str(date.today()) + '\n')
            fh.write('##Original File: ' + reference_annotation_basename + '\n')
            entries_written = 0
            for e in entries:
                try:
                    contig = str(e.get('contig', '.'))
                    # prefer reference coords if present
                    ref_pos = e.get('ref_pos', '')
                    add_pos = e.get('add_pos', '')
                    if ref_pos:
                        start, stop = ref_pos.split(',')
                        ftype = e.get('ref_type', '') or 'CDS'
                        source = reference_annotation_basename.split('_')[0] or 'reference'
                        info_attr = e.get('ref_info', '')
                    else:
                        # No ref pos, use add_pos coords
                        start, stop = add_pos.split(',') if add_pos else ('0','0')
                        ftype = e.get('add_type', '') or 'CDS'
                        source = e.get('add_type', '') or 'additional'
                        info_attr = e.get('add_info', '')
                    # attributes
                    attrs = []
                    attrs.append('Status=' + str(e.get('status', '')))
                    attrs.append('Coverage=' + str(e.get('coverage', '')))
                    if e.get('ref_info', ''):
                        attrs.append('Ref_info=' + str(e.get('ref_info', '')).replace(';','%3B'))
                    if e.get('add_info', ''):
                        attrs.append('Add_info=' + str(e.get('add_info', '')).replace(';','%3B'))
                    attr_str = ';'.join(attrs)
                    # construct GFF line
                    line = f"{contig}\t{source}\t{ftype}\t{start}\t{stop}\t.\t.\t.\t{attr_str}\n"
                    fh.write(line)
                    entries_written += 1
                except Exception:
                    # skip bad entry
                    continue
        logging.info('Wrote %d discordance GFF entries to %s', entries_written, report_path)
    except OSError as e:
        logging.error('Cannot write discordance GFF %s: %s', report_path, e)
        sys.exit(1)


def compute_discordance(ref_map_by_contig, add_map_by_contig, options):
    """Compare reference and additional maps per contig and return three lists:
    - only_in_ref: reference entries with no overlapping additional ORF
    - only_in_additional: additional ORFs that don't overlap any reference entry
    - mismatches: reference entries with overlapping additional ORFs that don't meet match criteria

    This version is strand-aware and will classify mismatches that are due to strand
    differences separately from type/coverage differences.

    Expected layouts:
    - ref entry: [strand, 'ref', type, info]
    - add entry: [strand, ..., type (index 3), info (last element)]
    """
    only_in_ref = []
    only_in_additional = []
    mismatches = []

    all_contigs = list(OrderedDict.fromkeys(list(ref_map_by_contig.keys()) + list(add_map_by_contig.keys())))
    matched_adds = set()
    cov_thresh = getattr(options, 'coverage', 100.0)

    for contig in all_contigs:
        ref_map = ref_map_by_contig.get(contig, {}) or {}
        add_map = add_map_by_contig.get(contig, {}) or {}

        # For each reference feature, find best overlapping additional ORF and classify
        for rpos, rdata in ref_map.items():
            rstart, rstop = _parse_pos(rpos)
            if rstart is None:
                continue
            rlen = (rstop - rstart + 1) if (rstop and rstart) else 0
            best_cov = 0.0
            best_add = None
            best_add_data = None
            matched = False

            # reference fields
            r_strand = rdata[0] if len(rdata) > 0 else ''
            r_type = rdata[3] if len(rdata) > 2 else ''
            r_info = rdata[-1] if len(rdata) > 3 else ''

            for apos, adata in add_map.items():
                astart, astop = _parse_pos(apos)
                if astart is None:
                    continue
                ov = max(0, min(rstop, astop) - max(rstart, astart) + 1)
                if ov <= 0:
                    continue
                cov = 100.0 * float(ov) / float(rlen) if rlen > 0 else 0.0
                if cov > best_cov:
                    best_cov = cov
                    best_add = apos
                    best_add_data = adata

                # additional fields
                a_strand = adata[0] if len(adata) > 0 else ''
                a_type = adata[3] if len(adata) > 3 else ''
                # frame check (distance of stops mod 3)
                try:
                    frame_ok = ((abs(astop - rstop) % 3) == 0)
                except Exception:
                    frame_ok = True

                # check for a fully satisfactory match: type, coverage, strand and frame
                if a_type == r_type and cov >= cov_thresh and (a_strand == r_strand) and frame_ok:
                    matched = True
                    matched_adds.add((contig, apos))
                    break

            if matched:
                # good match -> not discordant
                continue

            if best_add is None:
                # no overlapping additional ORF found
                only_in_ref.append({
                    'contig': contig,
                    'ref_pos': rpos,
                    'add_pos': '',
                    'ref_type': r_type,
                    'add_type': '',
                    'status': 'only_in_ref',
                    'coverage': '0.00',
                    'ref_info': r_info,
                    'add_info': ''
                })
            else:
                # overlapping additional ORF(s) exist but none satisfied the match criteria
                a_type = best_add_data[3] if len(best_add_data) > 3 else ''
                a_info = best_add_data[-1] if len(best_add_data) > 0 else ''
                a_strand = best_add_data[0] if len(best_add_data) > 0 else ''

                # compute reason flags
                type_match = (a_type == r_type)
                strand_match = (a_strand == r_strand)
                cov_ok = (best_cov >= cov_thresh)
                try:
                    # use frame between best add and ref
                    astart, astop = _parse_pos(best_add)
                    frame_ok = ((abs(astop - rstop) % 3) == 0) if (astop is not None) else True
                except Exception:
                    frame_ok = True

                # classify mismatch with strand-awareness
                if not cov_ok:
                    status = 'found_in_additional_but_below_coverage'
                elif not type_match and not strand_match:
                    status = 'found_in_additional_different_type_and_strand'
                elif not type_match:
                    status = 'found_in_additional_different_type'
                elif not strand_match:
                    status = 'found_in_additional_different_strand'
                elif not frame_ok:
                    status = 'found_in_additional_different_frame'
                else:
                    status = 'partial_overlap'

                mismatches.append({
                    'contig': contig,
                    'ref_pos': rpos,
                    'add_pos': best_add or '',
                    'ref_type': r_type,
                    'add_type': a_type,
                    'status': status,
                    'coverage': f"{best_cov:.2f}",
                    'ref_info': r_info,
                    'add_info': a_info,
                })

                if best_add:
                    matched_adds.add((contig, best_add))

        # Additional-only ORFs: those not matched and not overlapping any reference
        for apos, adata in add_map.items():
            if (contig, apos) in matched_adds:
                continue
            astart, astop = _parse_pos(apos)
            if astart is None:
                continue
            overlapped = False
            for rpos in ref_map.keys():
                rstart, rstop = _parse_pos(rpos)
                if rstart is None:
                    continue
                if max(rstart, astart) <= min(rstop, astop):
                    overlapped = True
                    break
            if not overlapped:
                only_in_additional.append({
                    'contig': contig,
                    'ref_pos': '',
                    'add_pos': apos,
                    'ref_type': '',
                    'add_type': adata[3] if len(adata) > 3 else '',
                    'status': 'only_in_additional',
                    'coverage': '0.00',
                    'ref_info': '',
                    'add_info': adata[-1] if len(adata) > 0 else '',
                })

    return only_in_ref, only_in_additional, mismatches


def comparator(options):
    # Multi-contig aware comparator
    genome_seq = ''
    genome_ID = None
    dna_regions = {}

    # Support both 'genome_DNA' and 'genome_dna' option names (compat with Annotation_Compare)
    genome_path = _get_opt(options, 'genome_DNA', 'genome_dna')

    # Load genome fasta if provided
    if genome_path:
        if not os.path.exists(genome_path):
            logging.error('Genome DNA file does not exist: %s', genome_path)
            sys.exit(1)
        try:
            fasta_in = gzip.open(genome_path, 'rt')
            dna_regions = fasta_load(fasta_in)
        except Exception:
            fasta_in = open(genome_path, 'r', encoding='unicode_escape')
            dna_regions = fasta_load(fasta_in)
        # genome_ID fallback
        try:
            if isinstance(dna_regions, dict) and len(dna_regions) > 0:
                genome_ID = next(iter(dna_regions.keys()))
                genome_seq = dna_regions[genome_ID]
            else:
                genome_ID = os.path.splitext(os.path.basename(genome_path))[0]
        except Exception:
            genome_seq = ''
            genome_ID = os.path.splitext(os.path.basename(genome_path))[0]
    else:
        # derive genome_ID from reference annotation filename
        genome_seq = ''
        genome_ID = os.path.splitext(os.path.basename(options.reference_annotation))[0]


    # Load reference annotation. If a tool-specific parser is requested, use it to ensure contig keys exist in dna_regions
    if getattr(options, 'reference_tool', None):
        try:
            reference_tool_mod = import_module('Tools.' + options.reference_tool + '.' + options.reference_tool,
                                               package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                reference_tool_mod = import_module(
                    'ORForise.Tools.' + options.reference_tool + '.' + options.reference_tool,
                    package='my_current_pkg')
            except ModuleNotFoundError:
                logging.error('Reference tool module not available: %s', options.reference_tool)
                sys.exit(1)
        reference_tool_fn = getattr(reference_tool_mod, options.reference_tool)
        try:
            # Call the tool parser; many tools return a contig->dict mapping. Ensure dna_regions contains those contig keys.
            ref_out = reference_tool_fn(options.reference_annotation, dna_regions)
            if isinstance(ref_out, dict):
                for contig_key in ref_out.keys():
                    if contig_key not in dna_regions:
                        dna_regions[contig_key] = ['']
        except Exception as e:
            logging.error('Failed to load reference annotation with tool %s: %s', options.reference_tool, e)
            sys.exit(1)
    else:
        try:
            gff_in = gzip.open(options.reference_annotation, 'rt')
            dna_regions = gff_load(options, gff_in, dna_regions)
        except Exception:
            gff_in = open(options.reference_annotation, 'r', encoding='unicode_escape')
            dna_regions = gff_load(options, gff_in, dna_regions)

    # Build ref_genes_by_contig: mapping contig -> OrderedDict(pos -> [strand, 'ref', type, info])
    ref_genes_by_contig = OrderedDict()

    if not getattr(options, 'reference_tool', None):
        # Parse GFF and group by seqid (first column)
        with open(options.reference_annotation, 'r', encoding='unicode_escape') as genome_gff:
            for line in genome_gff:
                line = line.rstrip('\n')
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                seqid = parts[0]
                ftype = parts[2]
                try:
                    # Determine if this feature type is requested
                    gene_types = options.gene_ident.split(',') if options.gene_ident else ['CDS']
                except Exception:
                    gene_types = ['CDS']
                if ftype not in gene_types and not ('CDS' in gene_types and ftype == 'CDS'):
                    # If user specified CDS and this is CDS, include; else skip
                    if ftype not in gene_types:
                        continue
                try:
                    start = int(parts[3])
                    stop = int(parts[4])
                    strand = parts[6]
                    pos = f"{start},{stop}"
                    info = parts[8]
                except (IndexError, ValueError):
                    continue
                if seqid not in ref_genes_by_contig:
                    ref_genes_by_contig[seqid] = OrderedDict()
                ref_genes_by_contig[seqid].update({pos: [strand, 'ref', ftype, info]})
    else:
        # Use a tool parser to produce ref_genes; expect tool to return mapping contig->dict
        try:
            reference_tool_mod = import_module('Tools.' + options.reference_tool + '.' + options.reference_tool,
                                               package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                reference_tool_mod = import_module('ORForise.Tools.' + options.reference_tool + '.' + options.reference_tool,
                                                   package='my_current_pkg')
            except ModuleNotFoundError:
                sys.exit("Tool not available")
        reference_tool_fn = getattr(reference_tool_mod, options.reference_tool)
        ref_out = reference_tool_fn(options.reference_annotation, dna_regions)
        # If the tool returns a mapping contig->dict, use that; otherwise assume single-contig and wrap
        if isinstance(ref_out, dict) and any(isinstance(v, dict) for v in ref_out.values()):
            ref_genes_by_contig = ref_out
        else:
            # single-contig output: place under genome_ID or first contig in dna_regions
            contig_key = genome_ID if genome_ID in dna_regions else (next(iter(dna_regions.keys())) if dna_regions else genome_ID)
            ref_genes_by_contig[contig_key] = ref_out

    # Get additional ORFs using tool parser; expect mapping contig->dict
    try:
        additional_tool_mod = import_module('Tools.' + options.additional_tool + '.' + options.additional_tool,
                                            package='my_current_pkg')
    except ModuleNotFoundError:
        try:
            additional_tool_mod = import_module('ORForise.Tools.' + options.additional_tool + '.' + options.additional_tool,
                                                package='my_current_pkg')
        except ModuleNotFoundError:
            sys.exit("Tool not available")
    additional_tool_fn = getattr(additional_tool_mod, options.additional_tool)
    additional_orfs = additional_tool_fn(options.additional_annotation, dna_regions)

    # Normalise additional_orfs: if single-contig dict, wrap under appropriate contig key
    if isinstance(additional_orfs, dict) and any(isinstance(v, dict) for v in additional_orfs.values()):
        additional_by_contig = additional_orfs
    else:
        contig_key = genome_ID if genome_ID in dna_regions else (next(iter(dna_regions.keys())) if dna_regions else genome_ID)
        additional_by_contig = {contig_key: additional_orfs}

    genes_To_Keep_by_contig = OrderedDict()

    # Iterate per contig and perform intersection logic
    for contig, orfs in additional_by_contig.items():
        ref_genes = ref_genes_by_contig.get(contig, OrderedDict())
        kept = OrderedDict()
        if options.coverage == 100.00:
            for orf, data in orfs.items():
                try:
                    o_Start = int(orf.split(',')[0])
                    o_Stop = int(orf.split(',')[1])
                except Exception:
                    continue
                o_Strand = data[0]
                additional_type = data[3]
                additional_info = data[-1]

                # Lookup exact-match reference entry safely
                ref_entry = ref_genes.get(f"{o_Start},{o_Stop}")
                if not ref_entry:
                    continue
                # ref_entry layout: [strand, 'ref', type, info]
                ref_type = ref_entry[3] if len(ref_entry) > 2 else ''
                ref_info = ref_entry[-1] if len(ref_entry) > 3 else ''

                if additional_type == ref_type and o_Strand == ref_entry[0]:
                    kept.update({f"{o_Start},{o_Stop}": [o_Strand, options.coverage, additional_type, ref_type, additional_info, ref_info]})
        else:
            cov_thresh = options.coverage
            for orf, data in orfs.items():
                try:
                    o_Start = int(orf.split(',')[0])
                    o_Stop = int(orf.split(',')[1])
                except Exception:
                    continue
                o_Strand = data[0]
                additional_type = data[3]
                additional_info = data[-1]

                for gene, r_data in ref_genes.items():
                    try:
                        g_Start = int(gene.split(',')[0])
                        g_Stop = int(gene.split(',')[1])
                    except Exception:
                        continue

                    # skip genes that start after this ORF (ref genes assumed sorted by start)
                    if g_Start > o_Stop:
                        break
                    # skip genes that end before this ORF
                    if g_Stop < o_Start:
                        continue

                    # compute overlap length without creating large sets
                    overlap = max(0, min(o_Stop, g_Stop) - max(o_Start, g_Start) + 1)
                    gene_len = (g_Stop - g_Start + 1)
                    if gene_len <= 0:
                        continue
                    cov = 100.0 * overlap / gene_len

                    g_Strand = r_data[0]
                    # r_data layout: [strand, 'ref', type, info]
                    ref_type = r_data[3] if len(r_data) > 2 else ''
                    ref_info = r_data[-1] if len(r_data) > 3 else ''

                    if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and cov >= cov_thresh:
                        if additional_type == ref_type:
                            kept[f"{g_Start},{g_Stop}"] = [g_Strand, int(cov), additional_type, ref_type,
                                                           additional_info, ref_info]
        genes_To_Keep_by_contig[contig] = sortORFs(kept)

    # Log counts for debugging why GFF might be empty
    try:
        total_ref = sum(len(v) for v in ref_genes_by_contig.values()) if ref_genes_by_contig else 0
    except Exception:
        total_ref = 0
    try:
        total_add = sum(len(v) for v in additional_by_contig.values()) if additional_by_contig else 0
    except Exception:
        total_add = 0
    try:
        total_kept = sum(len(v) for v in genes_To_Keep_by_contig.values()) if genes_To_Keep_by_contig else 0
    except Exception:
        total_kept = 0
    logging.info('Reference genes loaded: %d', total_ref)
    logging.info('Additional ORFs loaded: %d', total_add)
    logging.info('Kept genes after intersection: %d', total_kept)

    # If requested, compute discordance lists and write three GFF outputs
    if getattr(options, 'report_discordance', False):
        # Compute discordance lists
        only_in_ref, only_in_additional, mismatches = compute_discordance(ref_genes_by_contig, additional_by_contig, options)
        base = os.path.splitext(os.path.basename(options.output_file))[0] if getattr(options, 'output_file', None) else 'discordance'
        outdir = os.path.dirname(options.output_file) if getattr(options, 'output_file', None) else '.'
        ref_base = os.path.splitext(os.path.basename(options.reference_annotation))[0]

        # Keep the three detailed GFF outputs (backward compatible)
        gff_ref = os.path.join(outdir, f"{base}.only_in_reference.gff")
        gff_add = os.path.join(outdir, f"{base}.only_in_additional.gff")
        gff_mis = os.path.join(outdir, f"{base}.mismatches.gff")
        try:
            _write_discordance_gff(gff_ref, only_in_ref, ref_base)
            logging.info('Wrote discordance GFF: %s', gff_ref)
        except Exception:
            logging.exception('Failed to write discordance GFF: %s', gff_ref)
        try:
            _write_discordance_gff(gff_add, only_in_additional, ref_base)
            logging.info('Wrote discordance GFF: %s', gff_add)
        except Exception:
            logging.exception('Failed to write discordance GFF: %s', gff_add)
        try:
            _write_discordance_gff(gff_mis, mismatches, ref_base)
            logging.info('Wrote discordance GFF: %s', gff_mis)
        except Exception:
            logging.exception('Failed to write discordance GFF: %s', gff_mis)

        # Write a single concise summary TSV aggregating all discordance entries
        combined = []
        combined.extend(only_in_ref or [])
        combined.extend(mismatches or [])
        combined.extend(only_in_additional or [])
        combined_tsv = os.path.join(outdir, f"{base}.discordance_summary.tsv")
        try:
            _write_discordance_report(combined_tsv, combined)
            logging.info('Wrote discordance summary TSV: %s', combined_tsv)
        except Exception:
            logging.exception('Failed to write discordance summary TSV: %s', combined_tsv)

    # Ensure we always write a GFF (header + possibly diagnostic) so the core file is not empty
    genome_DNA_path = genome_path if genome_path else None



    # Write the kept genes GFF (this was missing and is why gff_writer wasn't called)
    try:
        logging.info('About to call gff_writer: total_kept=%d', total_kept)
        try:
            contig_summary = {c: len(v) for c, v in genes_To_Keep_by_contig.items()}
        except Exception:
            contig_summary = {}
        logging.info('Kept genes by contig (sample): %s', dict(list(contig_summary.items())[:10]))
        logging.info('Writing combined GFF to %s', options.output_file)
        gff_writer(genome_ID, genome_DNA_path, options.reference_annotation, getattr(options, 'reference_tool', None), None, options.additional_annotation, options.additional_tool, genes_To_Keep_by_contig, options.output_file, getattr(options, 'gene_ident', None))
        logging.info('gff_writer finished (check output file)')
    except Exception as e:
        logging.exception('Failed to write combined GFF: %s', e)

    # End of comparator


def main():
    print(WELCOME)

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': Annotation-Intersector Run Parameters')

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-ref', dest='reference_annotation', required=True,
                          help='Reference annotation GFF file')
    required.add_argument('-at', dest='additional_tool', required=True,
                          help='Tool name/format for additional annotation (module under Tools/)')
    required.add_argument('-add', dest='additional_annotation', required=True,
                          help='Additional annotation file to compare')
    required.add_argument('-o', dest='output_file', required=True,
                          help='Output GFF filename for kept genes')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-dna', dest='genome_DNA', required=False,
                          help='Genome DNA file (.fa) which both annotations are based on')
    optional.add_argument('-rt', dest='reference_tool', required=False,
                          help='Reference tool parser name (if not provided, GFF is expected)')
    optional.add_argument('-gi', dest='gene_ident', default='CDS', required=False,
                          help='Comma-separated feature types to consider from reference (default: CDS)')
    optional.add_argument('-cov', '--coverage', dest='coverage', default=100.0, type=float, required=False,
                          help='Percentage coverage threshold for intersection (default 100)')
    optional.add_argument('--report-discordance', dest='report_discordance', action='store_true', required=False,
                          help='If set, produce discordance reports (three GFFs)')
    optional.add_argument('--report-discordance-file', dest='report_discordance_file', required=False,
                          help='Optional base path for discordance reports')

    options = parser.parse_args()
    comparator(options)


if __name__ == '__main__':
    main()
    print('Complete')

