from svshare.svs import sv_end


# Insertions, breakends, and SVs without an end sit at their position.
def test_sv_end():
    assert sv_end("DEL", 100, 250) == 250
    assert sv_end("INS", 100, 5000) == 100
    assert sv_end("BND", 100, 9000) == 100
    assert sv_end("DUP", 100, None) == 100
    assert sv_end("INV", 100, 90) == 100
