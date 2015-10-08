# RVO: Rapid volume optimisation with an auxiliary equation of state #

[![DOI](https://zenodo.org/badge/14127/WMD-Bath/rvo.svg)](https://zenodo.org/badge/latestdoi/14127/WMD-Bath/rvo)

## Introduction ##

This is supporting information and a reference implemention of the RVO method.

## Contents ##

* **rvo.py** The main program file. May be used interactively; run `rvo.py -h` for usage information.

* **rvo_stats.py** Simulates the use of RVO with large sets of energy-volume curves.
Used to generate the figures for a research paper. Data files will be released on publication.
Run `rvo_stats.py -h` for usage information.

## Contents ##

* **dvxc.py** The main program file. May be used interactively; run `dvxc -h` for usage information.

* **dvxc_stats.py** Simulates the use of ΔV<sub>xc</sub> with large sets of energy-volume curves.
Used to generate the figures for an in-progress research paper. Data files will be released on publication. Run `dvxc -h` for usage information.

* **ase_generate_input.py** Generates an energy-volume curve for use with rvo.py, given a set of crystal structure files. Supports [file formats known to ASE](https://wiki.fysik.dtu.dk/ase/ase/io.html#module-ase.io).

* **interpolate_cell.py** Collection of functions used to interpolate lattice vectors to a target volume.

* **data/(files).dat** Cu4SnS4.dvxc.dat is provided as an example input file. This is an unpublished E-V curve for an interesting ternary chalcogenide.

* **interpolation_scheme.pdf** Supporting information; mathematical derivation of interpolation scheme used by **interpolate_cell.py**. The LaTeX source file is also provided.

## Requirements ##

Python 2.7 with Numpy, Scipy and Matplotlib.
Testing has largely been on unix-like filesystems; there may be some issues on Windows relating to file paths.

The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/) (ASE) is required for the "ase_generate_input" program.

## License ##

RVO is developed by the Walsh Materials Design group of the Department of Chemistry at the University of Bath, UK. Python code is licensed under the GNU General Public License (GPL) v3.
