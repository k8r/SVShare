# Read SV callers' VCFs and Jasmine's merge of them into memory.
from .svs import SV, Call, SVPosition, reference_end


# Yield (id, chrom, pos, info, format_values) for each record in a VCF, where
# info and format_values are dicts of the INFO and first sample's FORMAT fields.
def read_vcf_records(vcf):
    with open(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(
                item.split("=", 1) if "=" in item else (item, True)
                for item in fields[7].split(";")
            )
            fmt = {}
            if len(fields) > 9:
                fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
            yield fields[2], fields[0], int(fields[1]), info, fmt


def int_or_none(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


# For each record in Jasmine's merged VCF, yield (chrom, pos, vcf_end, svtype,
# svlen, members), where members lists the (input index, call ID) pairs merged
# into it. SUPP_VEC has one digit per input VCF, in file_list order, and IDLIST
# names the merged calls in that same order.
def read_jasmine_groups(merged_vcf):
    for _id, chrom, pos, info, _fmt in read_vcf_records(merged_vcf):
        inputs = [i for i, bit in enumerate(info["SUPP_VEC"]) if bit == "1"]
        members = list(zip(inputs, info["IDLIST"].split(",")))
        yield (chrom, pos, int_or_none(info.get("END")), info.get("SVTYPE"),
               int_or_none(info.get("SVLEN")), members)


# Return the calls in the given VCF, from the given sample and caller, by their
# VCF ID.
def read_calls(vcf, sample, caller):
    calls = {}
    for vcf_id, chrom, pos, info, fmt in read_vcf_records(vcf):
        calls[vcf_id] = Call(
            sample=sample, caller=caller, genotype=fmt.get("GT"),
            # Both sniffles2 and cuteSV report variant-supporting reads as DV.
            read_support=int_or_none(fmt.get("DV")), vcf_id=vcf_id, chrom=chrom, pos=pos,
            vcf_end=int_or_none(info.get("END")), svtype=info.get("SVTYPE"),
            svlen=int_or_none(info.get("SVLEN")),
        )
    return calls


# Return the SVs in the given Jasmine merged_vcf, each with its calls. The given
# calls_by_input is a list of dicts of calls by ID, one per input VCF in
# Jasmine's order.
def read_merged_svs(merged_vcf, calls_by_input):
    svs = []
    for sv_id, (chrom, pos, vcf_end, svtype, svlen, members) in enumerate(
        read_jasmine_groups(merged_vcf), start=1
    ):
        calls = []
        for index, vcf_id in members:
            calls.append(calls_by_input[index][vcf_id])
        position = SVPosition(
            chrom=chrom, pos=pos, vcf_end=vcf_end,
            reference_end=reference_end(svtype, pos, vcf_end), svtype=svtype, svlen=svlen,
        )
        svs.append(SV(sv_id, position, calls))
    return svs
