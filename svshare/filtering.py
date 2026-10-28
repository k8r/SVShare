# Functions for choosing which SVs a report includes.
import sqlite3

from .callers import CALLERS
from .population import DATASETS


# Combine --genes values and --gene-file lines into a sorted, deduplicated list.
def parse_gene_list(genes, gene_file):
    gene_set = set(genes or [])
    if gene_file:
        with open(gene_file) as f:
            gene_set.update(line.strip() for line in f if line.strip())
    return sorted(gene_set)


# The --caller-support and --min-samples filters, written in SQL. An SV passes
# if it meets caller_support in at least min_samples samples. Naming one caller
# also keeps SVs other callers called.
def _sample_support_condition(caller_support, min_samples):
    if caller_support == "any":
        samples = "SELECT sv_id, sample FROM calls GROUP BY sv_id, sample"
        params = []
    elif caller_support == "both":
        samples = ("SELECT sv_id, sample FROM calls GROUP BY sv_id, sample"
                   " HAVING COUNT(DISTINCT caller) = ?")
        params = [len(CALLERS)]
    else:
        samples = "SELECT sv_id, sample FROM calls WHERE caller = ? GROUP BY sv_id, sample"
        params = [caller_support]
    return (
        f"sv_id IN (SELECT sv_id FROM ({samples}) GROUP BY sv_id HAVING COUNT(*) >= ?)",
        params + [min_samples],
    )


# The --sv-type filter, written in SQL.
def _sv_type_condition(sv_type):
    if sv_type == "any":
        return None
    return "svtype = ?", [sv_type]


# The --max-frequency filter, written in SQL. An SV passes if, in each
# population database, its frequency is at most max_frequency or it isn't found there.
def _max_frequency_condition(max_frequency):
    if max_frequency is None:
        return None
    return (
        " AND ".join(f"({dataset}_af IS NULL OR {dataset}_af <= ?)" for dataset in DATASETS),
        [max_frequency] * len(DATASETS),
    )


# Return the IDs, in order, of the SVs in the database that pass the filters.
def filter_sv_ids(db_path, caller_support="any", min_samples=1, sv_type="any", max_frequency=None):
    conditions, params = [], []
    for condition in (
        _sample_support_condition(caller_support, min_samples),
        _sv_type_condition(sv_type),
        _max_frequency_condition(max_frequency),
    ):
        if condition:
            conditions.append(condition[0])
            params += condition[1]
    where = f" WHERE {' AND '.join(conditions)}" if conditions else ""
    with sqlite3.connect(db_path) as db:
        sv_ids = [sv_id for (sv_id,) in db.execute(f"SELECT sv_id FROM svs{where} ORDER BY sv_id", params)]
    db.close()
    return sv_ids
