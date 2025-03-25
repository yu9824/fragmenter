import math
import re
from collections.abc import Sequence

from rdkit import Chem

from fragmenter.typing import SMARTS


def calculate_fragment_scheme_order(
    fragmentation_scheme: dict[int, Sequence[SMARTS]],
) -> list[int]:
    """calculate `fragment_scheme_order` from `fragmentation_scheme`.

    Parameters
    ----------
    fragmentation_scheme : dict[int, Sequence[SMARTS]]
        fragmentation scheme

    Returns
    -------
    list[int]
        fragment_scheme_order
    """
    return list(
        map(
            lambda x: x[0],
            sorted(
                fragmentation_scheme.items(), key=_calculate_fragment_features
            )[::-1],
        )
    )


def _calculate_fragment_features(
    _tup_scheme: tuple[int, Sequence[SMARTS]],
) -> tuple[int, ...]:
    list_features: list[tuple[int, ...]] = list()
    for _smarts in _tup_scheme[1]:
        _mol_smarts = Chem.MolFromSmarts(_smarts)

        n_bonds = _mol_smarts.GetNumBonds()

        n_bonds_double = 0
        n_bonds_triple = 0
        for _bond in _mol_smarts.GetBonds():
            assert isinstance(_bond, Chem.Bond)
            _bond_type_dim = _bond.GetBondTypeAsDouble()
            n_bonds_double += math.isclose(_bond_type_dim, 2.0)
            n_bonds_triple += math.isclose(_bond_type_dim, 3.0)

        n_heavy_atoms = _mol_smarts.GetNumAtoms()
        n_total_atoms = n_heavy_atoms
        n_atoms_in_rings = 0
        n_atoms_aromatic = 0
        n_atoms_not_c_h = 0
        for _atom in _mol_smarts.GetAtoms():
            assert isinstance(_atom, Chem.QueryAtom)
            query = _atom.DescribeQuery()
            if _result_hydrogen_count := re.search(
                r"AtomHCount\s*(\d)", query
            ):
                n_total_atoms += int(_result_hydrogen_count.group(1))

            n_atoms_not_c_h += _atom.GetSymbol() not in {"C", "H"}
            n_atoms_aromatic += _atom.GetIsAromatic()

            n_atoms_in_rings += "AtomInNRings" in query

        list_features.append(
            (
                n_bonds,
                n_total_atoms,
                n_heavy_atoms,
                n_atoms_aromatic,
                n_atoms_not_c_h,
                n_bonds_double,
                n_bonds_triple,
            )
        )
    return max(list_features)


if __name__ == "__main__":
    from fragmenter.data import SMARTS_UNIFAC

    print(
        tuple(SMARTS_UNIFAC.items())[
            calculate_fragment_scheme_order(
                dict(enumerate(SMARTS_UNIFAC.values()))
            )[0]
        ]
    )
