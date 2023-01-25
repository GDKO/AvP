![Logo](assets/avp_logo.png)

# AvP

## Introduction

AvP performs automatic detection of HGT candidates within a phylogenetic framework. In a nutshell, the pipeline selects proteins based on a similarity hits file (blast or diamond), extracts sequences from the database, performs MSA, and phylogenetic inference. Then, it traverses the tree in order to detect candidate HGTs. Furthermore, it can evaluate the validity of candidate HGTs based on alternative topology testing and surrounding genes on the genome.


Documentation: [github wiki](https://github.com/GDKO/AvP/wiki)

## License and citation

The code is currently licensed under the GNU General Public License v3.0.

When using AvP, please cite [this paper](https://doi.org/10.1371/journal.pcbi.1010686):

Koutsovoulos GD, Granjeon Noriot S, Bailly-Bechet M, Danchin EGJ, Rancurel C (2022) 
**AvP: A software package for automatic phylogenetic detection of candidate horizontal gene transfers.**
*PLOS Computational Biology 18(11): e1010686*

doi:[10.1371/journal.pcbi.1010686](https://doi.org/10.1371/journal.pcbi.1010686)
