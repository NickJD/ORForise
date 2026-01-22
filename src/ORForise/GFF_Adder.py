from importlib import import_module
import argparse
from collections import OrderedDict, defaultdict, Counter
import gzip
from datetime import date
import sys

try:
    from .utils import *
except (ImportError, ModuleNotFoundError):
    from utils import *


########################################
def gff_writer(options, combined_ORFs_by_contig, output_file, reference_annotation, additional_annotation):
    write_out = open(output_file, 'w')
    write_out.write("##gff-version\t3\n#\tGFF-Adder\n#\tRun Date:" + str(date.today()) + '\n')
    write_out.write("##Genome DNA File:" + (options.genome_DNA if hasattr(options, 'genome_DNA') else '') + '\n')
    write_out.write("##Original File: " + reference_annotation + "\n##Additional File: " + additional_annotation + '\n')
    # meta counts
    Ref_Only = 0
    Ref_Combined = defaultdict(int)
    Non_Ref_Combined = defaultdict(int)

    # New counters: per-tool totals and per-contig per-tool breakdown
    per_tool_total = Counter()                       # counts for additional-only entries
    per_contig_per_tool = defaultdict(Counter)      # contig -> Counter(tool -> count)
    ref_per_tool_total = Counter()                  # counts for reference entries matched by tools
    ref_per_contig_per_tool = defaultdict(Counter)  # contig -> Counter(tool -> count)

    # Iterate contigs in deterministic order
    for contig in combined_ORFs_by_contig:
        combined_ORFs = combined_ORFs_by_contig[contig]
        # ref_gene_set for this contig: use keys from ref portion (we can detect by data[1]=='ref')
        ref_gene_set = [k for k, v in combined_ORFs.items() if len(v) > 1 and v[1] == 'ref']

        for pos, data in combined_ORFs.items():
            pos_ = pos.split(',')
            # pos may be like 'start,stop' or 'contig,start,stop' but here we expect 'start,stop'
            if len(pos_) >= 2:
                start = pos_[0]
                stop = pos_[-1]
            else:
                # fallback: skip malformed
                continue
            strand = data[0]

            # Build additional_annotation_info from the combined entry's additional list if present.
            # Normalise entries and prefer the info portion after any 'ToolName:info' prefix.
            additional_annotation_info = ''
            additional_items = []
            if len(data) > 4 and data[4]:
                for add in data[4]:
                    s = str(add).strip()
                    if not s:
                        continue
                    if ':' in s:
                        # split tool:info -> take info part
                        _, info_part = s.split(':', 1)
                        info_part = info_part.strip()
                    else:
                        info_part = s
                    if info_part:
                        additional_items.append(info_part)
                # deduplicate while preserving order
                seen = set()
                deduped = []
                for it in additional_items:
                    if it not in seen:
                        seen.add(it)
                        deduped.append(it)
                if deduped:
                    additional_annotation_info = ';'.join(deduped)
            elif len(data) > 3 and data[3]:
                additional_annotation_info = str(data[3]).strip()

            # tools list from options (maybe empty)
            tools = options.additional_tool.split(',') if getattr(options, 'additional_tool', None) else []

            # Determine matched_tools_list reliably:
            # prefer extracting tool names from data[4] entries, otherwise fallback to scanning values
            matched_tools_list = []
            try:
                if len(data) > 4 and data[4]:
                    for add in data[4]:
                        # expected format: "ToolName:info" or "info"
                        if isinstance(add, str) and ':' in add:
                            t = add.split(':', 1)[0].strip()
                            if t:
                                matched_tools_list.append(t)
                        else:
                            # try to detect one of the known tool names in the string
                            for tool in tools:
                                if tool and tool in str(add):
                                    matched_tools_list.append(tool)
                else:
                    # fallback: scan the whole data structure for tool names (previous behaviour)
                    for tool in tools:
                        if tool and any(tool in str(x) for x in data):
                            matched_tools_list.append(tool)
            except Exception:
                # keep matched_tools_list empty if anything unexpected happens
                matched_tools_list = []

            # Normalise matched tools: unique, sorted
            matched_tools_list = sorted(set(matched_tools_list))
            matched_tools_str = ','.join(matched_tools_list)
            matched_count = len(matched_tools_list)

            # Build GFF entry and update meta counters
            if pos not in ref_gene_set:  # Additional-only entry (not in reference)
                type_field = matched_tools_str if matched_tools_str else ''
                Non_Ref_Combined[matched_count] += 1

                # Update per-tool counters for additional-only entries
                if matched_tools_list:
                    for t in matched_tools_list:
                        per_tool_total[t] += 1
                        per_contig_per_tool[contig][t] += 1
                else:
                    # track unassigned additional entries (no tool name found)
                    per_tool_total['unassigned'] += 1
                    per_contig_per_tool[contig]['unassigned'] += 1

                if not getattr(options, 'clean', False):
                    entry = (contig + '\t' + type_field + '\t' + (data[2] if len(data) > 3 else '.') + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Additional_Annotations;' + (additional_annotation_info if additional_annotation_info else '') + '\n')
                else:
                    entry = (contig + '\t' + type_field + '\t' + (data[2] if len(data) > 3 else '.') + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\t' + (additional_annotation_info if additional_annotation_info else '') + '\n')

            else:
                # Reference entry
                if len(data) > 3 and data[3]:
                    info_field = data[3].replace('\n', '').strip()
                else:
                    info_field = '.'

                # Determine type and source fields
                type_field = data[2] if len(data) > 2 and data[2] else (
                    reference_annotation.split('/')[-1].split('.')[0] if reference_annotation else '.')
                source_field = data[5] if len(data) > 5 and data[5] else (
                    reference_annotation.split('/')[-1].split('.')[0] if reference_annotation else '')

                # If additional_annotation_info duplicates info_field content, remove duplicate fragments.
                filtered_additional = ''
                if additional_annotation_info:
                    add_parts = [p.strip() for p in str(additional_annotation_info).split(';') if p.strip()]
                    info_parts = [p.strip() for p in str(info_field).split(';') if p.strip() and p.strip() != '.']
                    filtered = []
                    for ap in add_parts:
                        dup = False
                        for ip in info_parts:
                            # treat duplication if exact match or obvious substring relationship
                            if ip and (ap == ip or ip in ap or ap in ip):
                                dup = True
                                break
                        if not dup:
                            filtered.append(ap)
                    filtered_additional = ';'.join(filtered)

                if not filtered_additional:
                    # Reference-only (no meaningful unique additional annotations)
                    Ref_Only += 1
                    if not getattr(options, 'clean', False):
                        entry = (
                                    contig + '\t' + source_field + '\t' + type_field + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation;' + info_field + '\n')
                    else:
                        entry = (
                                    contig + '\t' + source_field + '\t' + type_field + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\t' + info_field + '\n')
                else:
                    # Reference entry that had additional annotations (combined)
                    Ref_Combined[matched_count] += 1

                    # Update per-tool counters for reference-matched entries
                    if matched_tools_list:
                        for t in matched_tools_list:
                            ref_per_tool_total[t] += 1
                            ref_per_contig_per_tool[contig][t] += 1
                    else:
                        ref_per_tool_total['unassigned'] += 1
                        ref_per_contig_per_tool[contig]['unassigned'] += 1

                    if not getattr(options, 'clean', False):
                        entry = (
                                    contig + '\t' + source_field + '\t' + type_field + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation;' + info_field + ';Matched_Annotations=' + filtered_additional + '\n')
                    else:
                        entry = (
                                    contig + '\t' + source_field + '\t' + type_field + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\t' + info_field + '\n')

            write_out.write(entry)

    # Produce metadata output if requested
    if getattr(options, 'output_meta', False) == True:
        # Summaries
        total_ref_combined = sum(Ref_Combined.values())
        total_nonref = sum(Non_Ref_Combined.values())
        total_reference_genes = Ref_Only + total_ref_combined
        Ref_Combined_counter = Counter(Ref_Combined)
        Non_Ref_Combined_counter = Counter(Non_Ref_Combined)
        with open(output_file.replace('.gff','_Meta.txt'),'w') as meta_out:
            meta_out.write("GFF-Adder Metadata Report\n")
            meta_out.write("=========================\n")
            meta_out.write("Run Date: {}\n".format(date.today()))
            meta_out.write("Genome DNA: {}\n".format(getattr(options, 'genome_DNA', 'N/A')))
            meta_out.write("Reference annotation: {}\n".format(reference_annotation))
            meta_out.write("Additional annotation(s): {}\n\n".format(additional_annotation))

            meta_out.write("Summary counts\n")
            meta_out.write("--------------\n")
            meta_out.write(f"Reference-only genes (no matching additional annotation): {Ref_Only}\n")
            meta_out.write(f"Reference genes with additional matches (combined): {total_ref_combined}\n")
            meta_out.write(f"TOTAL reference genes observed: {total_reference_genes}\n")
            meta_out.write(f"Additional-only genes (not present in reference): {total_nonref}\n\n")

            meta_out.write("Distribution of matches for reference-combined entries (num matched tools -> count)\n")
            if Ref_Combined_counter:
                for k in sorted(Ref_Combined_counter):
                    meta_out.write(f"  {k:>3}    {Ref_Combined_counter[k]}\n")
            else:
                meta_out.write("  None\n")
            meta_out.write("\n")

            meta_out.write("Distribution of matches for non-reference entries (num matched tools -> count)\n")
            if Non_Ref_Combined_counter:
                for k in sorted(Non_Ref_Combined_counter):
                    meta_out.write(f"  {k:>3}    {Non_Ref_Combined_counter[k]}\n")
            else:
                meta_out.write("  None\n")
            meta_out.write("\n")

            # Per-tool totals for additional-only entries
            meta_out.write("Per-tool additional-only annotation totals\n")
            meta_out.write("-----------------------------------------\n")
            if per_tool_total:
                for t, c in per_tool_total.most_common():
                    meta_out.write(f"  {t}: {c}\n")
            else:
                meta_out.write("  None\n")
            meta_out.write("\n")

            # Per-contig breakdown for additional-only entries: only contigs with additional genes
            meta_out.write("Per-contig breakdown for additional-only annotations (only contigs with additions shown)\n")
            meta_out.write("---------------------------------------------------------------------------------------\n")
            if per_contig_per_tool:
                for contig in sorted(per_contig_per_tool):
                    counter = per_contig_per_tool[contig]
                    if sum(counter.values()) == 0:
                        continue
                    meta_out.write(f"  {contig}:\n")
                    for t, c in counter.most_common():
                        meta_out.write(f"    {t}: {c}\n")
                    meta_out.write("\n")
            else:
                meta_out.write("  None\n\n")

            # Per-tool totals and per-contig breakdown for reference genes matched by additional tools
            meta_out.write("Per-tool totals for reference genes matched by additional tools\n")
            meta_out.write("---------------------------------------------------------------\n")
            if ref_per_tool_total:
                for t, c in ref_per_tool_total.most_common():
                    meta_out.write(f"  {t}: {c}\n")
            else:
                meta_out.write("  None\n")
            meta_out.write("\n")

            meta_out.write("Per-contig breakdown for reference genes matched by additional tools\n")
            meta_out.write("--------------------------------------------------------------------\n")
            if ref_per_contig_per_tool:
                for contig in sorted(ref_per_contig_per_tool):
                    counter = ref_per_contig_per_tool[contig]
                    if sum(counter.values()) == 0:
                        continue
                    meta_out.write(f"  {contig}:\n")
                    for t, c in counter.most_common():
                        meta_out.write(f"    {t}: {c}\n")
                    meta_out.write("\n")
            else:
                meta_out.write("  None\n\n")

            meta_out.write("Notes\n")
            meta_out.write("-----\n")
            meta_out.write(" - 'Reference-only' are reference entries that had no recorded additional annotation information.\n")
            meta_out.write(" - 'Per-tool' counts are based on the tool names extracted from the additional-annotation provenance\n")
            meta_out.write("   (expected format in combined entries: 'ToolName:info'). Entries with no tool detected are shown as 'unassigned'.\n")
            meta_out.write("\nEnd of report\n")



def gff_adder(options):
    # Load fasta into dna_regions (supports multi-contig)
    try:
        try:
            fasta_in = gzip.open(options.genome_DNA, 'rt')
            dna_regions = fasta_load(fasta_in)
        except Exception:
            fasta_in = open(options.genome_DNA, 'r', encoding='unicode_escape')
            dna_regions = fasta_load(fasta_in)
    except Exception:
        # Fallback to legacy single-contig behaviour
        genome_seq = ""
        with open(options.genome_DNA, 'r') as genome_fasta:
            for line in genome_fasta:
                line = line.replace("\n", "")
                if not line.startswith('>'):
                    genome_seq += str(line)
                else:
                    genome_ID = line.split()[0].replace('>','')
        # Create dna_regions with single entry
        dna_regions = OrderedDict()
        dna_regions[genome_ID] = (genome_seq, len(genome_seq), list(), None)

    ###########################################
    # Build reference gene dict per-contig
    ref_genes_by_contig = defaultdict(OrderedDict)

    if not options.reference_tool:  # IF using Ensembl/file for comparison
        # Parse reference gff to populate ref_genes_by_contig (retain original info fields)
        # Detect gzip by magic bytes (first two bytes)
        is_gz = False
        with open(options.reference_annotation, 'rb') as _probe:
            magic = _probe.read(2)
            is_gz = (magic == b'\x1f\x8b')

        if is_gz:
            gff_in = gzip.open(options.reference_annotation, 'rt', errors='replace')
        else:
            # Open as plain text, replace undecodable bytes rather than fail
            gff_in = open(options.reference_annotation, 'r', encoding='utf-8', errors='replace')

        count = 0
        try:
            for line in gff_in:
                if line.startswith('#') or line.strip() == '':
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                contig = parts[0]
                if contig not in dna_regions:
                    # skip records for contigs not in provided fasta
                    continue
                try:
                    if 'CDS' in options.gene_ident and len(options.gene_ident) == 1:
                        if "CDS" in parts[2] and len(parts) == 9:
                            start = int(parts[3])
                            stop = int(parts[4])
                            strand = parts[6]
                            pos = f"{start},{stop}"
                            # store as [strand, source, type, info] to match downstream expectations
                            info = parts[8]
                            ref_genes_by_contig[contig][pos] = [strand, parts[1], parts[2], info]
                            count += 1
                    else:
                        gene_types = options.gene_ident.split(',')
                        if any(gene_type in parts[2] for gene_type in gene_types):
                            start = int(parts[3])
                            stop = int(parts[4])
                            strand = parts[6]
                            pos = f"{start},{stop}"
                            # store as [strand, source, type, info] to match downstream expectations
                            info = parts[8]
                            ref_genes_by_contig[contig][pos] = [strand, parts[1], parts[2], info]
                            count += 1
                except IndexError:
                    continue
        finally:
            try:
                gff_in.close()
            except Exception:
                pass

    else:
        # Reference tool provided: attempt to call it with dna_regions first (multi-contig aware), fallback to legacy signature
        reference_tool = options.reference_tool if options.reference_tool != 'StORF-Reporter' else 'StORF-Reporter'
        try:
            reference_tool_mod = import_module('Tools.' + reference_tool + '.' + reference_tool, package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                reference_tool_mod = import_module('ORForise.Tools.' + reference_tool + '.' + reference_tool, package='my_current_pkg')
            except ModuleNotFoundError:
                sys.exit("Tool not available")
        reference_tool_func = getattr(reference_tool_mod, reference_tool)
        # Try multi-contig signature
        try:
            ref_result = reference_tool_func(options.reference_annotation, dna_regions)
        except TypeError:
            # Fallback to legacy signature, try passing genome seq string
            genome_seq = ''.join([dna_regions[c][0] for c in dna_regions])
            ref_result = reference_tool_func(reference_annotation=options.reference_annotation, genome_seq=genome_seq, gene_ident=options.gene_ident)
        # Expect ref_result to be dict of contig -> {pos: data}
        for contig, mapping in ref_result.items() if isinstance(ref_result, dict) else []:
            ref_genes_by_contig[contig].update(mapping)

    # Ensure each contig has an OrderedDict even if empty
    for contig in dna_regions:
        if contig not in ref_genes_by_contig:
            ref_genes_by_contig[contig] = OrderedDict()

    # Collect additional annotations per contig
    additional_annotations_by_contig = defaultdict(OrderedDict)
    tool_count = 0
    for tool in options.additional_tool.split(','):
        try:
            additional_tool_mod = import_module('Tools.' + tool + '.' + tool, package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                additional_tool_mod = import_module('ORForise.Tools.' + tool + '.' + tool, package='my_current_pkg')
            except ModuleNotFoundError:
                sys.exit("Tool not available")
        additional_tool_func = getattr(additional_tool_mod, tool)

        anno_file = options.additional_annotation.split(',')[tool_count]
        tool_count += 1
        # Try calling tool in multi-contig mode first
        try:
            tool_orfs = additional_tool_func(anno_file, dna_regions)
        except TypeError:
            # Fallback to legacy signature expecting genome_seq
            genome_seq = ''.join([dna_regions[c][0] for c in dna_regions])
            tool_orfs = additional_tool_func(anno_file, genome_seq, options.gene_ident)

        # tool_orfs may be either {contig: {pos: data}} or a flat {pos: data}
        if isinstance(tool_orfs, dict):
            # If top-level keys look like contig names (present in dna_regions) then treat as multi-contig
            top_keys = list(tool_orfs.keys())
            if top_keys and top_keys[0] in dna_regions:
                for contig, mapping in tool_orfs.items():
                    # Merge mapping into additional_annotations_by_contig and record tool provenance
                    for pos_k, pos_v in mapping.items():
                        # store tuple (value, tool)
                        additional_annotations_by_contig[contig][pos_k] = (pos_v, tool)
            else:
                # Treat as flat mapping — assume single contig if only one contig present
                if len(dna_regions) == 1:
                    only_contig = next(iter(dna_regions))
                    for pos_k, pos_v in tool_orfs.items():
                        additional_annotations_by_contig[only_contig][pos_k] = (pos_v, tool)
                else:
                    # If multiple contigs but mapping has contig-prefixed keys like 'contig,start,stop', split them
                    for k, v in tool_orfs.items():
                        parts = k.split(',')
                        if len(parts) == 3 and parts[0] in dna_regions:
                            contig = parts[0]
                            pos = parts[1] + ',' + parts[2]
                            additional_annotations_by_contig[contig][pos] = (v, tool)
                        else:
                            # Unknown format: assign nowhere (skip)
                            continue
        tool_orfs = None

    # Combine per-contig: keep reference entries and append additional annotations as supplemental
    combined_ORFs_by_contig = OrderedDict()
    for contig in dna_regions:
        combined = OrderedDict()
        # Add reference entries first; normalise to [strand, 'ref', type, ref_info, additional_list]
        for pos, val in ref_genes_by_contig.get(contig, {}).items():
            strand = val[0] if len(val) > 0 else '.'
            src = 'ref'
            source_field = val[1] if len(val) > 1 else 'ref'
            ftype = val[2] if len(val) > 2 else '.'
            ref_info = val[3] if len(val) > 3 else '.'
            combined[pos] = [strand, src, ftype, ref_info, [], source_field]

        # Now incorporate additional annotations without overwriting reference entries
        for pos, wrapped in additional_annotations_by_contig.get(contig, {}).items():
            # wrapped is (value, tool)
            if isinstance(wrapped, tuple) and len(wrapped) == 2:
                val, toolname = wrapped
            else:
                val = wrapped
                toolname = ''

            # Extract strand/type/info from value heuristically
            strand_a = val[0] if isinstance(val, (list, tuple)) and len(val) > 0 else '.'
            ftype_a = val[3] if isinstance(val, (list, tuple)) and len(val) > 2 else '.'
            info_a = ''
            if isinstance(val, (list, tuple)) and len(val) > 3:
                info_a = val[4]
            elif isinstance(val, str):
                info_a = val

            # If matching pos exists in reference, append additional info to its additional list
            if pos in combined:
                addstr = (toolname + ':' + info_a) if toolname else info_a
                combined[pos][4].append(addstr)
            else:
                # Create a new entry for additional-only annotation: [strand, 'add', type, '.', [tool:info]]
                addstr = (toolname + ':' + info_a) if toolname else info_a
                #combined[pos] = [strand_a, 'add', ftype_a if ftype_a else '.', '.', [addstr]]
                combined[pos] = [strand_a, 'add', ftype_a if ftype_a else '.', '.', [addstr], toolname]
        # Sort ORFs for this contig
        combined = sortORFs(combined)
        combined_ORFs_by_contig[contig] = combined

    # Call writer
    gff_writer(options, combined_ORFs_by_contig, options.output_file, options.reference_annotation, options.additional_annotation)


def main():
    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': GFF-Adder Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-dna', dest='genome_DNA', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                    'are based on')
    required.add_argument('-ref', dest='reference_annotation', required=True,
                        help='Which reference annotation file to use as reference?')
    required.add_argument('-at', dest='additional_tool', required=True,
                        help='Which format to use for additional annotation? - Can provide multiple annotations (Tool1,Tool2)')
    required.add_argument('-add', dest='additional_annotation', required=True,
                        help='Which annotation file to add to reference annotation? - Can provide multiple annotations (1.GFF,2.GFF)')
    required.add_argument('-o', dest='output_file', required=True,
                        help='Output filename')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rt', dest='reference_tool', required=False,
                        help='Which tool format to use as reference? - If not provided, will default to the '
                             'standard GFF format and will only look for "CDS" features')
    optional.add_argument('--gene_ident', action='store', dest='gene_ident', default='CDS',
                        help='Identifier used for identifying genomic features in reference annotation "CDS,rRNA,tRNA"')
    optional.add_argument('-mc', dest='mark_consensus', action='store_true', required=False,
                        help='Default - False: Mark reference annotations which where present in the additional tool annotation')
    optional.add_argument('-c', dest='clean', action='store_true', required=False,
                        help='Default - False: Do not mark 9th column with "Original/Matched/Additional tag"')
    optional.add_argument('--meta', dest='output_meta', action='store_true', required=False,
                        help='Default - False: Output metadata file')
    optional.add_argument('--olap', dest='overlap', default=50, type=int, required=False,
                        help='Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50 nt')

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-v', dest='verbose', default='False', type=eval, choices=[True, False],
                      help='Default - False: Print out runtime status')

    options = parser.parse_args()

    gff_adder(options)



if __name__ == "__main__":
    try:
        try:
            main()
        except Exception:
            logging.exception('Unhandled exception in main')
    finally:
        print(CLOSING)