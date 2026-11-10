# Build small databases for tests, and read them back.
from svshare.database import SCHEMA, connect


# Write a database with the given rows, each a dict of column values, and
# return its path.
def make_database(tmp_path, svs=(), calls=(), transcripts=(), sv_transcript_parts=(), genes=()):
    db_path = tmp_path / "svshare.db"
    tables = {"svs": svs, "calls": calls, "transcripts": transcripts,
              "sv_transcript_parts": sv_transcript_parts, "genes": genes}
    with connect(db_path) as db:
        db.executescript(SCHEMA)
        for table, rows in tables.items():
            for row in rows:
                db.execute(f"INSERT INTO {table} ({', '.join(row)})"
                           f" VALUES ({', '.join('?' * len(row))})", list(row.values()))
    return db_path


# Return the given columns, or all of them, of each row in table as a tuple, in
# the order the rows were saved.
def read_rows(db_path, table, *columns):
    with connect(db_path) as db:
        return [tuple(row) for row in db.execute(
            f"SELECT {', '.join(columns) or '*'} FROM {table} ORDER BY rowid"
        )]
