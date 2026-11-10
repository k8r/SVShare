import sqlite3

import pytest

from svshare.database import connect


# Changes are saved, and the database is closed even when something fails.
def test_connect_saves_and_closes(tmp_path):
    db_path = tmp_path / "test.db"
    with connect(db_path) as db:
        db.execute("CREATE TABLE t (x INTEGER)")
        db.execute("INSERT INTO t VALUES (1)")
    with pytest.raises(sqlite3.ProgrammingError):
        db.execute("SELECT 1")

    with pytest.raises(ValueError):
        with connect(db_path) as db:
            db.execute("INSERT INTO t VALUES (2)")
            raise ValueError
    with pytest.raises(sqlite3.ProgrammingError):
        db.execute("SELECT 1")
    with connect(db_path) as db:
        assert db.execute("SELECT x FROM t").fetchall() == [(1,)]
