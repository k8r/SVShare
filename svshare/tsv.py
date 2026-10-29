# Write a report's SVs as a TSV table, one row per SV.
import csv

from .callers import CALLERS

COLUMNS = [
    "chrom", "pos", "end", "svtype", "svlen", "samples", "callers", "gnomad_af",
    "colorsdb_af", "genes", "gene_names", "haploinsufficiency", "triplosensitivity",
]


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


# Join each gene's value as "GENE: value", skipping genes without one.
def _per_gene(genes, field):
    return "; ".join(f"{gene.name}: {getattr(gene, field)}" for gene in genes
                     if getattr(gene, field)) or None


# Return an SV's COLUMNS values.
def _row(sv):
    return {
        "chrom": sv.chrom,
        "pos": sv.pos,
        "end": sv.end,
        "svtype": sv.svtype,
        "svlen": sv.svlen,
        "samples": _samples(sv.calls),
        "callers": ",".join(caller for caller in CALLERS
                            if any(call.caller == caller for call in sv.calls)),
        "gnomad_af": sv.gnomad_af,
        "colorsdb_af": sv.colorsdb_af,
        "genes": ", ".join(f"{gene.name} ({', '.join(gene.regions)})" for gene in sv.genes) or None,
        "gene_names": _per_gene(sv.genes, "full_name"),
        "haploinsufficiency": _per_gene(sv.genes, "haploinsufficiency"),
        "triplosensitivity": _per_gene(sv.genes, "triplosensitivity"),
    }


# Write SVs to output as a TSV table.
def write_tsv(svs, output):
    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(_row(sv) for sv in svs)
