# `report` subcommand: filter and report analysis results.
import sys

from .callers import CALLERS
from .database import SVFilter, select_sv_ids, sv_calls, sv_rows, sv_transcript_parts
from .population import DATASETS
from .svs import SV, Call, GeneDetails, SVPosition, Transcript, TranscriptPart
from .tsv import write_tsv

# Combine --genes values and --gene-file lines into a sorted, deduplicated list.
def parse_gene_list(genes, gene_file):
    gene_set = set(genes or [])
    if gene_file:
        with open(gene_file) as f:
            gene_set.update(line.strip() for line in f if line.strip())
    return sorted(gene_set)


# Turn the report's filter options into an SVFilter.
def make_sv_filter(caller_support, min_samples, sv_type, max_frequency, genes,
                   overlaps_any_gene, region, near_reference_gap):
    if caller_support == "any":
        callers = ()
    elif caller_support == "both":
        callers = tuple(CALLERS)
    else:
        callers = (caller_support,)
    return SVFilter(
        callers=callers,
        min_samples=min_samples,
        svtype=None if sv_type == "any" else sv_type,
        population_datasets=tuple(DATASETS),
        max_allele_frequency=max_frequency,
        genes=tuple(genes or ()),
        overlaps_any_transcript=overlaps_any_gene,
        region=region,
        near_reference_gap={"include": None, "only": True, "exclude": False}[near_reference_gap],
    )


# Sort key for chromosomes: numbered ones in number order, then X, Y, and M,
# then the rest by name.
def _chrom_order(chrom):
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name), "")
    if name in ("X", "Y", "M"):
        return (1, "XYM".index(name), "")
    return (2, 0, name)


# Return the SVs in sv_ids, in genome order.
def load_svs(db_path, sv_ids):
    calls = {}
    for row in sv_calls(db_path, sv_ids):
        calls.setdefault(row["sv_id"], []).append(
            Call(row["sample"], row["caller"], row["genotype"], row["read_support"],
                 row["reference_support"])
        )
    overlapped, nearby = _load_transcript_parts(db_path, sv_ids)

    svs = []
    for row in sv_rows(db_path, sv_ids):
        position = SVPosition(
            chrom=row["chrom"], pos=row["pos"], vcf_end=row["vcf_end"],
            reference_end=row["reference_end"], svtype=row["svtype"], svlen=row["svlen"],
        )
        svs.append(SV(
            row["sv_id"], position, calls.get(row["sv_id"], []),
            gnomad_af=row["gnomad_af"], colorsdb_af=row["colorsdb_af"],
            near_reference_gap=row["near_reference_gap"],
            overlapped_parts=overlapped.get(row["sv_id"], []),
            nearby_parts=nearby.get(row["sv_id"], []),
        ))
    svs.sort(key=lambda sv: (_chrom_order(sv.position.chrom), sv.position.pos, sv.sv_id))
    return svs


# Return the transcript parts each SV in sv_ids overlaps and the ones only near
# it, each by SV ID. Parts of the same transcript share one Transcript, and
# transcripts of the same gene share one GeneDetails.
def _load_transcript_parts(db_path, sv_ids):
    genes, transcripts, overlapped, nearby = {}, {}, {}, {}
    for row in sv_transcript_parts(db_path, sv_ids):
        gene_id, transcript_id = row["gene_id"], row["transcript_id"]
        if gene_id not in genes:
            genes[gene_id] = GeneDetails(
                gene_id, row["hgnc_id"], row["full_name"], row["summary"],
                row["haploinsufficiency_score"], row["haploinsufficiency"],
                row["triplosensitivity_score"], row["triplosensitivity"],
            )
        if transcript_id not in transcripts:
            transcripts[transcript_id] = Transcript(
                transcript_id, row["gene_name"], gene_id, row["transcript_type"],
                bool(row["mane_select"]), row["chrom"], row["strand"], genes[gene_id],
            )
        part = TranscriptPart(
            transcripts[transcript_id], row["region"], row["number"], row["pos"],
            row["end_pos"],
        )
        parts = overlapped if row["overlaps_sv"] else nearby
        parts.setdefault(row["sv_id"], []).append(part)
    return overlapped, nearby


def run(args):
    # SQLite would otherwise create an empty database at a mistyped path.
    if not args.database.is_file():
        sys.exit(f"svshare report: {args.database} isn't a file; pass the database "
                 "analyze wrote, e.g. results/svshare.db")
    genes = parse_gene_list(args.genes, args.gene_file)
    output = args.output or (args.database.parent / f"report.{args.format}")
    print(f"Reporting on {args.database}")
    print(f"  caller-support: {args.caller_support}")
    print(f"  sv-type: {args.sv_type}")
    print(f"  max-frequency: {args.max_frequency}")
    print(f"  min-samples: {args.min_samples}")
    print(f"  near-reference-gap: {args.near_reference_gap}")
    if genes:
        print(f"  genes: {', '.join(genes)}")
    if args.overlaps_any_gene:
        print("  overlaps-any-gene: yes")
    if args.region:
        chrom, start, end = args.region
        print(f"  region: {chrom}:{start}-{end}")

    sv_filter = make_sv_filter(caller_support=args.caller_support,
                               min_samples=args.min_samples, sv_type=args.sv_type,
                               max_frequency=args.max_frequency, genes=genes,
                               overlaps_any_gene=args.overlaps_any_gene,
                               region=args.region,
                               near_reference_gap=args.near_reference_gap)
    sv_ids = select_sv_ids(args.database, sv_filter)
    print(f"{len(sv_ids)} SVs pass the filters")
    svs = load_svs(args.database, sv_ids)
    if genes:
        missing = []
        for gene in genes:
            overlapped = False
            for sv in svs:
                if sv.overlaps_gene(gene):
                    overlapped = True
            if not overlapped:
                missing.append(gene)
        if missing:
            print(f"  No SVs that pass the filters overlap: {', '.join(missing)}")
    if args.format == "tsv":
        write_tsv(svs, output)
        print(f"Wrote {output}")
    else:
        print("HTML reports aren't written yet; pass --format tsv for a table")
