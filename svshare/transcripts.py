# Find the transcripts near each SV and the parts of them each SV affects.
import gzip
import re
import shlex
import sqlite3
from pathlib import Path

from .commands import run_command

# The GTF doesn't include promoters, so each transcript's promoter is taken to be
# the PROMOTER_SIZE bp just upstream of its start.
PROMOTER_SIZE = 2000

# How much of each transcript to keep around an SV, for showing it in context:
# parts within VIEW_MARGIN_FRACTION of the SV's length on either side, but at
# least MIN_VIEW_MARGIN bp.
VIEW_MARGIN_FRACTION = 0.25
MIN_VIEW_MARGIN = 1000

# GTF features always copied into BED rows.
CODING_FEATURES = {"CDS", "UTR"}

# Files written to the output directory.
SV_WINDOWS_BED = "sv_windows.bed"
TRANSCRIPT_REGIONS_BED = "transcript_regions.bed"
OVERLAPS_BED = "sv_transcript_overlaps.bed"

# Matches one GTF attribute, quoted (gene_name "BRCA1";) or not (exon_number 3;).
_GTF_ATTRIBUTE = re.compile(r'(\S+) "?([^";]*)"?;')


# Yield (chrom, start, end, feature, strand, attributes, tags) for each line of
# a basic transcript in a GENCODE GTF. Basic transcripts are the ones GENCODE
# tags as full-length and well supported. start and end are 1-based, as in the GTF.
def read_basic_records(gtf):
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as f:
        for line in f:
            # Checking the raw line first skips the rest of the file cheaply.
            if 'tag "basic"' not in line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            attributes, tags = {}, set()
            for key, value in _GTF_ATTRIBUTE.findall(fields[8]):
                if key == "tag":
                    tags.add(value)
                else:
                    attributes[key] = value
            yield (fields[0], int(fields[3]), int(fields[4]), fields[2], fields[6],
                   attributes, tags)


# Return BED rows (chrom, start, end, gene_name, gene_id, transcript_id,
# transcript_type, mane_select, region, number, strand) for the parts of each
# basic transcript in the GTF.
def transcript_regions(gtf):
    rows = []
    exons = {}
    codons = {}
    coding = set()
    for chrom, start, end, feature, strand, attributes, tags in read_basic_records(gtf):
        transcript_id = attributes["transcript_id"]
        transcript = (attributes["gene_name"], attributes["gene_id"], transcript_id,
                      attributes["transcript_type"], int("MANE_Select" in tags))
        if feature in CODING_FEATURES:
            coding.add(transcript_id)
            rows.append((chrom, start - 1, end, *transcript, feature,
                         int(attributes["exon_number"]), strand))
        elif feature == "exon":
            exons.setdefault(transcript_id, []).append(
                (chrom, start, end, int(attributes["exon_number"]), transcript, strand)
            )
        elif feature in ("start_codon", "stop_codon"):
            # Saved until every line is read, then cut out of the CDS or UTR
            # containing them.
            codons.setdefault(transcript_id, []).append(
                (chrom, start - 1, end, *transcript, feature, int(attributes["exon_number"]),
                 strand)
            )
        elif feature == "transcript":
            # The promoter is upstream of the transcript's start, which is its
            # right-hand end on the minus strand.
            if strand == "+":
                promoter = (max(0, start - 1 - PROMOTER_SIZE), start - 1)
            else:
                promoter = (end, end + PROMOTER_SIZE)
            rows.append((chrom, *promoter, *transcript, "promoter", None, strand))

    rows = split_out_codons(rows, codons)
    for transcript_id, transcript_exons in exons.items():
        transcript_exons.sort()
        # Non-coding transcripts have no CDS or UTR lines, so their exons are
        # stored instead.
        if transcript_id not in coding:
            rows.extend((chrom, start - 1, end, *transcript, "noncoding_exon", number, strand)
                        for chrom, start, end, number, transcript, strand in transcript_exons)
        # Introns are the gaps between exons, numbered like the exon before them.
        for left, right in zip(transcript_exons, transcript_exons[1:]):
            chrom, _, left_end, left_number, transcript, strand = left
            right_start, right_number = right[1], right[3]
            if right_start - 1 > left_end:
                rows.append((chrom, left_end, right_start - 1, *transcript, "intron",
                             min(left_number, right_number), strand))
    return rows


# Return rows with each codon cut out of the CDS or UTR row containing it and
# added as its own row, so its bases aren't labeled twice. codons holds each
# transcript's codon rows, keyed by transcript ID.
def split_out_codons(rows, codons):
    result = []
    for row in rows:
        cuts = codons.get(row[5], []) if row[8] in CODING_FEATURES else []
        start, end = row[1], row[2]
        for _, cut_start, cut_end, *_ in sorted(cuts, key=lambda cut: cut[1]):
            if cut_start < end and cut_end > start:
                if start < cut_start:
                    result.append((row[0], start, cut_start, *row[3:]))
                start = max(start, cut_end)
        if start < end:
            result.append((row[0], start, end, *row[3:]))
    for transcript_codons in codons.values():
        result.extend(transcript_codons)
    return result


# Write a BED file with each SV in the database widened by its view margin.
def write_sv_windows_bed(db, path):
    with open(path, "w") as f:
        for sv_id, chrom, pos, end_pos, svtype in db.execute(
            "SELECT sv_id, chrom, pos, end_pos, svtype FROM svs"
        ):
            start = pos - 1
            if svtype in ("INS", "BND") or end_pos is None:
                # An insertion or breakend covers just its position, whatever END
                # its caller gave, so its window is the minimum margin on each side.
                end = pos
            else:
                end = max(end_pos, pos)
            margin = max(MIN_VIEW_MARGIN, round(VIEW_MARGIN_FRACTION * (end - start)))
            window_start, window_end = max(0, start - margin), end + margin
            row = (chrom, window_start, window_end, sv_id, start, end)
            f.write("\t".join(map(str, row)) + "\n")


# Intersect the SV windows with the transcript regions from the GTF and save the
# overlaps to the database. regions_bed, if given, is reused, or built if it
# doesn't exist. Returns (bedtools command, SVs hitting a transcript, total SVs).
def add_sv_transcripts(db_path, gtf, out_dir, regions_bed=None):
    out_dir = Path(out_dir)
    windows_bed = out_dir / SV_WINDOWS_BED
    overlaps_bed = out_dir / OVERLAPS_BED
    reuse = regions_bed is not None and Path(regions_bed).exists()
    regions_bed = Path(regions_bed) if regions_bed else out_dir / TRANSCRIPT_REGIONS_BED

    if not reuse:
        with open(regions_bed, "w") as f:
            f.writelines(
                "\t".join("." if value is None else str(value) for value in row) + "\n"
                for row in transcript_regions(gtf)
            )

    with sqlite3.connect(db_path) as db:
        write_sv_windows_bed(db, windows_bed)
        cmd = ["bedtools", "intersect", "-a", str(windows_bed), "-b", str(regions_bed),
               "-wa", "-wb"]
        with open(overlaps_bed, "w") as f:
            run_command(cmd, stdout=f)
        command = f"{shlex.join(cmd)} > {shlex.quote(str(overlaps_bed))}"

        # Each line is the SV window's 6 columns followed by the transcript
        # region's 11. Parts are stored with 1-based, inclusive positions, like
        # the rest of the database.
        transcripts, parts = {}, set()
        with open(overlaps_bed) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                sv_id, sv_start, sv_end = int(fields[3]), int(fields[4]), int(fields[5])
                (chrom, start, end, gene_name, gene_id, transcript_id, transcript_type,
                 mane_select, region, number, strand) = fields[6:17]
                start, end = int(start), int(end)
                transcripts[transcript_id] = (gene_name, gene_id, transcript_type,
                                              int(mane_select), chrom, strand)
                parts.add((sv_id, transcript_id, region,
                           None if number == "." else int(number), start + 1, end,
                           int(start < sv_end and end > sv_start)))

        db.executemany(
            "INSERT INTO transcripts (transcript_id, gene_name, gene_id, transcript_type,"
            " mane_select, chrom, strand) VALUES (?, ?, ?, ?, ?, ?, ?)",
            [(transcript_id, *values) for transcript_id, values in sorted(transcripts.items())],
        )
        db.executemany(
            "INSERT INTO sv_transcript_parts (sv_id, transcript_id, region, number, pos,"
            " end_pos, overlaps_sv) VALUES (?, ?, ?, ?, ?, ?, ?)",
            sorted(parts, key=lambda part: (part[0], part[1], part[4])),
        )
        db.execute(
            "INSERT INTO commands (sample, tool, command) VALUES (NULL, 'bedtools', ?)",
            (command,),
        )
        total = db.execute("SELECT COUNT(*) FROM svs").fetchone()[0]
    db.close()
    affected = {part[0] for part in parts if part[6]}
    return command, len(affected), total
