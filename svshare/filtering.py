# Functions for choosing which SVs a report includes.
import sqlite3

from .callers import CALLERS


# Combine --genes values and --gene-file lines into a sorted, deduplicated list.
def parse_gene_list(genes, gene_file):
    gene_set = set(genes or [])
    if gene_file:
        with open(gene_file) as f:
            gene_set.update(line.strip() for line in f if line.strip())
    return sorted(gene_set)


# The --caller-support filter, written in SQL. A caller counts as calling an SV
# if it called it in any sample, and naming one caller also keeps SVs other
# callers called.
def _caller_support_condition(caller_support):
    if caller_support == "any":
        return None
    if caller_support == "both":
        return (
            "sv_id IN (SELECT sv_id FROM calls GROUP BY sv_id"
            " HAVING COUNT(DISTINCT caller) = ?)",
            [len(CALLERS)],
        )
    return "sv_id IN (SELECT sv_id FROM calls WHERE caller = ?)", [caller_support]


# The --sv-type filter, written in SQL.
def _sv_type_condition(sv_type):
    if sv_type == "any":
        return None
    return "svtype = ?", [sv_type]


# Return the IDs, in order, of the SVs in the database that pass the filters.
def filter_sv_ids(db_path, caller_support="any", sv_type="any"):
    conditions, params = [], []
    for condition in (_caller_support_condition(caller_support), _sv_type_condition(sv_type)):
        if condition:
            conditions.append(condition[0])
            params += condition[1]
    where = f" WHERE {' AND '.join(conditions)}" if conditions else ""
    with sqlite3.connect(db_path) as db:
        sv_ids = [sv_id for (sv_id,) in db.execute(f"SELECT sv_id FROM svs{where} ORDER BY sv_id", params)]
    db.close()
    return sv_ids
