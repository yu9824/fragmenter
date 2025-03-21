# FRAGMENTER

## Install

```bash
pip install git+https://github.com/yu9824/fragmenter.git
```

***This repository is forked from [simonmb/fragmentation_algorithm](https://github.com/simonmb/fragmentation_algorithm).***

> # fragmentation_algorithm (new version)
>
> ![Algorithm for automatically fragmenting molecules into predefined groups](imgs/graphical%20abstract2_white_bg.png)
>
> - This algorithm allows fragmenting molecules in an automatic manner based on groups defined as SMARTS.
> - A file is included which executes the fragmentation as was done in the paper and one that shows how to use it regularly.
> - Dependencies are rdkit, tqdm (optional for the fragmentation as was done in the paper)
> - Why is there a newer version?
>   - The fragmentation algorithm originally published with [the paper](https://doi.org/10.1186/s13321-019-0382-3) tried to find smaller functional groups that are contained within other larger functional groups to sort them automatically in an intelligent way. This turned out to be quite difficult to implement as the capabilities of RDKit to match SMARTS with SMARTS are  limited. This lead me to writing workarounds that became broken in subsequent RDKit versions.
>   - The newer version does not try to automatically sort the groups you are searching for but relies on the user to provide this order. This allows the algorithm to be applicable with recent RDKit versions without problems.
> - Added AI Docstrings and capability to draw molecules with groups to enhance debugging the creation of SMARTS, here is one example:
>
> <p align="center">
>   <img src="imgs/example2.png" alt="example of molecule with found groups" width="500"/>
> </p>
>
> ---
> The [original version of this algorithm](https://github.com/simonmb/fragmentation_algorithm_paper) was published with a paper applying it for UNIFAC:
> - Title: **Flexible Heuristic Algorithm for Automatic Molecule Fragmentation: Application to the UNIFAC Group Contribution Model**
> - DOI: [10.1186/s13321-019-0382-39](https://doi.org/10.1186/s13321-019-0382-3)
