# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: base311
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
from fragmenter_utils import draw_mol_with_highlights_and_legend, get_table_with_atom_properties_relevant_to_SMARTS
from rdkit import Chem
import SMARTS

UNIFAC_SMARTS = SMARTS.UNIFAC.copy()

fragmentation_scheme = {i+1: j[1] for i, j in enumerate(UNIFAC_SMARTS)}
simple_fragmenter = fragmenter(fragmentation_scheme, algorithm='simple')

group_names = {k+1:v[0] for k,v  in enumerate(UNIFAC_SMARTS)}
for i, SMILES in enumerate(['CC1=C(C(=CC(=C1Cl)Cl)Cl)Cl', 'CCC=CCC1=C(CCC1=O)Cc1ccccc1', 'CCC1=CC(=C(C=C1)CC)CC']):
        mol = Chem.MolFromSmiles(SMILES)
        fragmentation, success, fragmentation_matches = simple_fragmenter.fragment(mol)
        img = draw_mol_with_highlights_and_legend(mol, fragmentation_matches, group_names)
        _, _, _, formatted_rows = get_table_with_atom_properties_relevant_to_SMARTS(mol)
        display(img)
        print('\n'.join(formatted_rows))
        
