# Write a report's SVs as a TSV table, one row per SV.
import csv

from .callers import CALLERS

COLUMNS = [
    "genome_order", "genes", "overlap_summary", "svtype", "svlen", "gnomad_af",
    "colorsdb_af", "chrom", "pos", "end", "haploinsufficiency", "triplosensitivity",
    "samples", "callers", "weakest_sample_support", "read_support",
    "near_reference_gap", "gene_full_names",
]

# The order a gene's parts are listed in, roughly from upstream to downstream.
# Any other kind of part is listed after these, alphabetically.
REGION_ORDER = [
    "promoter", "UTR", "start_codon", "CDS", "stop_codon", "noncoding_transcript_exon",
    "intron",
]

# The kind of part each transcript region counts as in the overlap_summary column.
REGION_KINDS = {
    "CDS": "exon", "UTR": "exon", "start_codon": "exon", "stop_codon": "exon",
    "noncoding_transcript_exon": "exon", "intron": "intron", "promoter": "promoter",
}


def _region_order(region):
    if region in REGION_ORDER:
        return (REGION_ORDER.index(region), "")
    return (len(REGION_ORDER), region)


# List the samples an SV was called in, each with its genotype.
def _samples(calls):
    genotypes = {}
    for call in calls:
        genotypes.setdefault(call.sample, {})
        if call.genotype and "." not in call.genotype:
            genotypes[call.sample][call.caller] = call.genotype
    samples = []
    for sample, by_caller in sorted(genotypes.items()):
        if not by_caller:
            samples.append(sample)
        elif len(set(by_caller.values())) == 1:
            samples.append(f"{sample} ({next(iter(by_caller.values()))})")
        else:
            listed = ", ".join(f"{by_caller[caller]} {caller}" for caller in CALLERS
                               if caller in by_caller)
            samples.append(f"{sample} ({listed})")
    return ", ".join(samples)


# List how many reads supported each call, out of all the reads its caller
# counted there if it counted the reference ones, by sample and then caller,
# leaving out calls without a count.
def _read_support(calls):
    calls = sorted((call for call in calls if call.read_support is not None),
                   key=lambda call: (call.sample, list(CALLERS).index(call.caller)))
    supports = []
    for call in calls:
        support = str(call.read_support)
        if call.reference_support is not None:
            support += f" of {call.read_support + call.reference_support}"
        supports.append(f"{support} from {call.caller} for {call.sample}")
    return "; ".join(supports) or None


# The read support of the least supported sample, taking each sample's
# best-supported call.
def _weakest_sample_support(calls):
    best = {}
    for call in calls:
        if call.read_support is not None:
            best[call.sample] = max(best.get(call.sample, 0), call.read_support)
    return min(best.values(), default=None)


# Return (name, details, regions) for each gene an SV overlaps, by name, where
# regions are the kinds of part it overlaps across the gene's transcripts.
def _overlapped_genes(sv):
    genes = {}
    for part in sv.overlapped_parts:
        transcript = part.transcript
        _details, regions = genes.setdefault(
            (transcript.gene_name, transcript.gene_id), (transcript.gene_details, set())
        )
        regions.add(part.region)
    return [(name, details, sorted(regions, key=_region_order))
            for (name, _gene_id), (details, regions) in sorted(genes.items())]


# List the kinds of part an SV overlaps, merged across all its genes, in
# alphabetical order so the column sorts into groups. The genes column has each
# gene's parts.
def _overlap_summary(sv):
    kinds = set()
    for part in sv.overlapped_parts:
        kinds.add(REGION_KINDS.get(part.region, part.region))
    return ", ".join(sorted(kinds)) or None


# Join each gene's detail as "GENE: value", skipping genes without one.
def _per_gene(genes, field):
    values = [(name, getattr(details, field) if details else None) for name, details, _ in genes]
    return "; ".join(f"{name}: {value}" for name, value in values if value) or None


# Return an SV's COLUMNS values.
def _row(sv):
    genes = _overlapped_genes(sv)
    return {
        "chrom": sv.position.chrom,
        "pos": sv.position.pos,
        "end": sv.position.reference_end,
        "svtype": sv.position.svtype,
        "svlen": sv.position.svlen,
        "samples": _samples(sv.calls),
        "read_support": _read_support(sv.calls),
        "weakest_sample_support": _weakest_sample_support(sv.calls),
        "near_reference_gap": sv.near_reference_gap,
        "callers": ",".join(caller for caller in CALLERS
                            if any(call.caller == caller for call in sv.calls)),
        "gnomad_af": sv.gnomad_af,
        "colorsdb_af": sv.colorsdb_af,
        "genes": ", ".join(f"{name} ({', '.join(regions)})" for name, _, regions in genes)
                 or None,
        "overlap_summary": _overlap_summary(sv),
        "gene_full_names": _per_gene(genes, "full_name"),
        "haploinsufficiency": _per_gene(genes, "haploinsufficiency"),
        "triplosensitivity": _per_gene(genes, "triplosensitivity"),
    }


# Write SVs to output as a TSV table, numbering them in the order given so the
# table can be sorted back into that order.
def write_tsv(svs, output):
    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({"genome_order": i, **_row(sv)}
                         for i, sv in enumerate(svs, start=1))
