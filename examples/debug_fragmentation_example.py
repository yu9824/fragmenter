# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: py311
#     language: python
#     name: python3
# ---

# %% [markdown]
# Debugging how SMARTS should look like can be quite challenging. THis is why I developed these fragmeneter_utils. On another note, it always helps to have a look on how rdkit implements the properties, sometimes things are unfortunately not so easy to grasp. I have opened several issues thinking they were bugs that turned out to be hard chemically to solve or design decisions. Especially, I think aromaticity is hard. (At least for me)
# You can get more info from here:
# https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
# https://greglandrum.github.io/rdkit-blog/posts/2023-05-26-drawing-options-explained.html
#

# %%
from IPython.display import display

from fragmenter import fragmenter
from fragmenter.utils import (
    calculate_fragment_scheme_order,
    draw_mol_with_highlights_and_legend,
    get_table_with_atom_properties_relevant_to_SMARTS,
)
from rdkit import Chem
from fragmenter.data import SMARTS_UNIFAC

# %%
fragmentation_scheme = {
    index: list(tup_smarts)
    for index, tup_smarts in enumerate(SMARTS_UNIFAC.values(), start=1)
}
fragmentation_scheme_order = calculate_fragment_scheme_order(
    fragmentation_scheme
)
group_names = {
    index: name for index, name in enumerate(SMARTS_UNIFAC.keys(), start=1)
}
simple_fragmenter = fragmenter(
    fragmentation_scheme, fragmentation_scheme_order, algorithm="simple"
)
# simple_fragmenter = fragmenter(
#     fragmentation_scheme,
#     fragmentation_scheme_order,
#     algorithm="complete",
#     function_to_choose_fragmentation=lambda x: x[0],
# )

for i, smiles in enumerate(
    [
        "CC1=C(C(=CC(=C1Cl)Cl)Cl)Cl",
        "CCC=CCC1=C(CCC1=O)Cc1ccccc1",
        "CCC1=CC(=C(C=C1)CC)CC",
    ]
):
    mol = Chem.MolFromSmiles(smiles)
    fragmentation, success, fragmentation_matches = simple_fragmenter.fragment(
        mol
    )

    img = draw_mol_with_highlights_and_legend(
        mol, fragmentation_matches, group_names
    )
    _, _, _, formatted_rows = (
        get_table_with_atom_properties_relevant_to_SMARTS(mol)
    )
    display(img)
    print("\n".join(formatted_rows))

