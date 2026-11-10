# Build small databases for tests.
import sqlite3

from svshare.database import SCHEMA


# Write a database with the given rows, each a dict of column values, and
# return its path.
def make_database(tmp_path, svs=(), calls=(), transcripts=(), sv_transcript_parts=(), genes=()):
    db_path = tmp_path / "svshare.db"
    tables = {"svs": svs, "calls": calls, "transcripts": transcripts,
              "sv_transcript_parts": sv_transcript_parts, "genes": genes}
    with sqlite3.connect(db_path) as db:
        db.executescript(SCHEMA)
        for table, rows in tables.items():
            for row in rows:
                db.execute(f"INSERT INTO {table} ({', '.join(row)})"
                           f" VALUES ({', '.join('?' * len(row))})", list(row.values()))
    db.close()
    return db_path
