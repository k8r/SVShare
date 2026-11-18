import gzip

from svshare.svs import SV, SVPosition, Transcript, TranscriptPart
from svshare.transcripts import (
    count_affected_svs, find_sv_transcripts, transcript_regions, write_sv_windows_bed,
)

MANE = 'tag "basic"; tag "MANE_Select";'
BASIC = 'tag "basic";'

# A GENCODE-style GTF, with 1-based, inclusive coordinates, and four transcripts:
#   T1   GENE1's MANE Select transcript, on the + strand, with two exons, UTRs at
#        both ends and an intron between the exons
#   T1b  a non-coding GENE1 transcript with one exon inside T1's intron
#   T2   GENE2's MANE Select transcript, on the - strand, with two coding exons
#        numbered from the right, since that's where it starts
#   T1x  a GENE1 transcript that isn't basic, which should be ignored
# As in GENCODE, each start codon is inside a CDS and T1's stop codon is inside
# its UTR. The last field is the exon number, which transcript lines don't have.
GTF_RECORDS = [
    ("transcript", 1001, 5000, "+", "GENE1", "T1", "protein_coding", MANE, None),
    ("exon", 1001, 2000, "+", "GENE1", "T1", "protein_coding", MANE, 1),
    ("exon", 3001, 5000, "+", "GENE1", "T1", "protein_coding", MANE, 2),
    ("UTR", 1001, 1200, "+", "GENE1", "T1", "protein_coding", MANE, 1),
    ("CDS", 1201, 2000, "+", "GENE1", "T1", "protein_coding", MANE, 1),
    ("start_codon", 1201, 1203, "+", "GENE1", "T1", "protein_coding", MANE, 1),
    ("CDS", 3001, 4800, "+", "GENE1", "T1", "protein_coding", MANE, 2),
    ("UTR", 4801, 5000, "+", "GENE1", "T1", "protein_coding", MANE, 2),
    ("stop_codon", 4801, 4803, "+", "GENE1", "T1", "protein_coding", MANE, 2),
    ("transcript", 2201, 2400, "+", "GENE1", "T1b", "lncRNA", BASIC, None),
    ("exon", 2201, 2400, "+", "GENE1", "T1b", "lncRNA", BASIC, 1),
    ("transcript", 2201, 2300, "+", "GENE1", "T1x", "retained_intron", "", None),
    ("exon", 2201, 2300, "+", "GENE1", "T1x", "retained_intron", "", 1),
    ("transcript", 20001, 21000, "-", "GENE2", "T2", "protein_coding", MANE, None),
    ("exon", 20601, 21000, "-", "GENE2", "T2", "protein_coding", MANE, 1),
    ("exon", 20001, 20400, "-", "GENE2", "T2", "protein_coding", MANE, 2),
    ("CDS", 20601, 21000, "-", "GENE2", "T2", "protein_coding", MANE, 1),
    ("start_codon", 20998, 21000, "-", "GENE2", "T2", "protein_coding", MANE, 1),
    ("CDS", 20001, 20400, "-", "GENE2", "T2", "protein_coding", MANE, 2),
]

# The fields each transcript adds to its BED rows: gene_name, gene_id,
# transcript_id, transcript_type, and mane_select.
T1 = ("GENE1", "GENE1_ID", "T1", "protein_coding", 1)
T1B = ("GENE1", "GENE1_ID", "T1b", "lncRNA", 0)
T2 = ("GENE2", "GENE2_ID", "T2", "protein_coding", 1)
T3 = ("GENE3", "GENE3_ID", "T3", "protein_coding", 1)


def write_gtf(tmp_path, records=GTF_RECORDS):
    path = tmp_path / "genes.gtf.gz"
    with gzip.open(path, "wt") as f:
        f.write("##description: test annotation\n")
        for feature, start, end, strand, gene, transcript, type_, tags, exon in records:
            exon_number = "" if exon is None else f"exon_number {exon}; "
            f.write(
                f"chr1\tHAVANA\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
                f'gene_id "{gene}_ID"; transcript_id "{transcript}"; gene_name "{gene}"; '
                f'transcript_type "{type_}"; {exon_number}{tags}\n'
            )
    return path


# Checks every BED row built from GTF_RECORDS. The cases it covers:
#   - codons are cut out of the CDS or UTR containing them: T1's start and stop
#     codons, and T2's start codon, at the right end of its CDS on the - strand
#   - introns take the number of the exon before them in reading order, so T2's
#     is 1, since its exon 1 is on the right
#   - T1b, which is non-coding, gets an exon row instead of CDS and UTR rows
#   - T1x, which isn't basic, gets no rows
#   - promoters are the 2 kb upstream of each transcript: left of T1 and T1b,
#     with T1's cut off at the chromosome start, and right of T2
def test_transcript_regions_covers_each_part_of_basic_transcripts(tmp_path):
    assert sorted(transcript_regions(write_gtf(tmp_path)), key=str) == sorted([
        ("chr1", 0, 1000, *T1, "promoter", None, "+"),
        ("chr1", 1000, 1200, *T1, "UTR", 1, "+"),
        ("chr1", 1200, 1203, *T1, "start_codon", 1, "+"),
        ("chr1", 1203, 2000, *T1, "CDS", 1, "+"),
        ("chr1", 2000, 3000, *T1, "intron", 1, "+"),
        ("chr1", 3000, 4800, *T1, "CDS", 2, "+"),
        ("chr1", 4800, 4803, *T1, "stop_codon", 2, "+"),
        ("chr1", 4803, 5000, *T1, "UTR", 2, "+"),
        ("chr1", 200, 2200, *T1B, "promoter", None, "+"),
        ("chr1", 2200, 2400, *T1B, "noncoding_exon", 1, "+"),
        ("chr1", 20000, 20400, *T2, "CDS", 2, "-"),
        ("chr1", 20400, 20600, *T2, "intron", 1, "-"),
        ("chr1", 20600, 20997, *T2, "CDS", 1, "-"),
        ("chr1", 20997, 21000, *T2, "start_codon", 1, "-"),
        ("chr1", 21000, 23000, *T2, "promoter", None, "-"),
    ], key=str)


# A GTF with one transcript, T3, on the + strand, whose stop codon is split by
# its intron:
#   exon 1  CDS, then the stop codon's first two bases, which GENCODE labels UTR
#   exon 2  the stop codon's last base, then the rest of the UTR
SPLIT_STOP_RECORDS = [
    ("transcript", 1001, 4000, "+", "GENE3", "T3", "protein_coding", MANE, None),
    ("exon", 1001, 2000, "+", "GENE3", "T3", "protein_coding", MANE, 1),
    ("exon", 3001, 4000, "+", "GENE3", "T3", "protein_coding", MANE, 2),
    ("CDS", 1001, 1998, "+", "GENE3", "T3", "protein_coding", MANE, 1),
    ("UTR", 1999, 2000, "+", "GENE3", "T3", "protein_coding", MANE, 1),
    ("UTR", 3001, 4000, "+", "GENE3", "T3", "protein_coding", MANE, 2),
    ("stop_codon", 1999, 2000, "+", "GENE3", "T3", "protein_coding", MANE, 1),
    ("stop_codon", 3001, 3001, "+", "GENE3", "T3", "protein_coding", MANE, 2),
]


# Checks every BED row built from SPLIT_STOP_RECORDS. The cases it covers:
#   - each piece of the stop codon is cut out of its own exon's UTR
#   - exon 1's UTR, which is all stop codon, leaves no row behind
#   - exon 2's UTR is shortened by one base
def test_transcript_regions_splits_a_codon_split_across_exons(tmp_path):
    assert sorted(transcript_regions(write_gtf(tmp_path, SPLIT_STOP_RECORDS)),
                  key=str) == sorted([
        ("chr1", 0, 1000, *T3, "promoter", None, "+"),
        ("chr1", 1000, 1998, *T3, "CDS", 1, "+"),
        ("chr1", 1998, 2000, *T3, "stop_codon", 1, "+"),
        ("chr1", 2000, 3000, *T3, "intron", 1, "+"),
        ("chr1", 3000, 3001, *T3, "stop_codon", 2, "+"),
        ("chr1", 3001, 4000, *T3, "UTR", 2, "+"),
    ], key=str)


# Checks which transcript parts two deletions get, one at 1501-2500 and one away
# from all genes. The cases it covers:
#   - the first deletion's window is 501-3500, the 1 kb minimum margin on each
#     side, so it gets every part of T1 and T1b in it, and nothing of T2
#   - parts have their numbers and 1-based positions, ordered by transcript and
#     then position
#   - the parts the deletion itself overlaps, T1's first CDS piece, T1's intron,
#     T1b's promoter, and T1b's exon, are kept apart from the ones only nearby
#   - the second deletion gets no parts
def test_find_sv_transcripts_finds_parts_around_each_sv(tmp_path):
    svs = [
        SV(1, SVPosition("chr1", 1501, 2500, 2500, "DEL", -1000), []),
        SV(2, SVPosition("chr1", 50001, 50100, 50100, "DEL", -100), []),
    ]

    find_sv_transcripts(svs, write_gtf(tmp_path), tmp_path)

    t1 = Transcript("T1", "GENE1", "GENE1_ID", "protein_coding", True, "chr1", "+")
    t1b = Transcript("T1b", "GENE1", "GENE1_ID", "lncRNA", False, "chr1", "+")
    assert svs[0].overlapped_parts == [
        TranscriptPart(t1, "CDS", 1, 1204, 2000),
        TranscriptPart(t1, "intron", 1, 2001, 3000),
        TranscriptPart(t1b, "promoter", None, 201, 2200),
        TranscriptPart(t1b, "noncoding_exon", 1, 2201, 2400),
    ]
    assert svs[0].nearby_parts == [
        TranscriptPart(t1, "promoter", None, 1, 1000),
        TranscriptPart(t1, "UTR", 1, 1001, 1200),
        TranscriptPart(t1, "start_codon", 1, 1201, 1203),
        TranscriptPart(t1, "CDS", 2, 3001, 4800),
    ]
    assert svs[1].overlapped_parts == []
    assert svs[1].nearby_parts == []
    assert count_affected_svs(svs) == 1


# Checks that an existing regions BED is reused instead of rebuilt from the GTF.
def test_find_sv_transcripts_reuses_existing_regions_bed(tmp_path):
    svs = [SV(1, SVPosition("chr1", 1501, 1600, 1600, "DEL", -100), [])]
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("\t".join(map(str, ("chr1", 1200, 2000, *T1, "CDS", 1, "+"))) + "\n")

    find_sv_transcripts(svs, tmp_path / "missing.gtf.gz", tmp_path, regions_bed)

    t1 = Transcript("T1", "GENE1", "GENE1_ID", "protein_coding", True, "chr1", "+")
    assert svs[0].overlapped_parts == [TranscriptPart(t1, "CDS", 1, 1201, 2000)]


# Checks that each SV's window is built around where it ends on the reference,
# not its VCF END.
def test_write_sv_windows_bed_uses_reference_end(tmp_path):
    svs = [
        SV(1, SVPosition("chr1", 5000, 5001, 5000, "INS", 300), []),
        SV(2, SVPosition("chr1", 9001, 13000, 13000, "DEL", -4000), []),
    ]
    windows_bed = tmp_path / "windows.bed"

    write_sv_windows_bed(svs, windows_bed)

    assert windows_bed.read_text().splitlines() == [
        "chr1\t3999\t6000\t1\t4999\t5000",
        "chr1\t8000\t14000\t2\t9000\t13000",
    ]
