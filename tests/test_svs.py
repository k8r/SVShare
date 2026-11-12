from svshare.svs import reference_end


# Insertions, breakends, and SVs without an end sit at their position.
def test_reference_end():
    assert reference_end("DEL", 100, 250) == 250
    assert reference_end("INS", 100, 5000) == 100
    assert reference_end("BND", 100, 9000) == 100
    assert reference_end("DUP", 100, None) == 100
    assert reference_end("INV", 100, 90) == 100
