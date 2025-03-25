from fragmenter import fragmenter
from fragmenter.data import SMARTS_UNIFAC

fragmentation_scheme = dict(enumerate(SMARTS_UNIFAC.values()))
fragmentation_scheme_order = dict(enumerate(SMARTS_UNIFAC.keys()))


def test_simple():
    fragmenter_simple = fragmenter(
        fragmentation_scheme, fragmentation_scheme_order, algorithm="simple"
    )
    fragmenter_simple.fragment("C1CCCCC1")


def test_complete():
    fragmenter_simple = fragmenter(
        fragmentation_scheme,
        fragmentation_scheme_order,
        algorithm="complete",
        n_heavy_atoms_cuttoff=10,
    )
    fragmenter_simple.fragment("C1CCCCC1")
