# README

## Overview
This repo gives a solution to processesing crystal structure data and creates molecular subgraphs based on atomic connectivity. It extracts the ligands from the crystal structure and visualizes both the molecule, line graphs and 2nd-order linegraph with detailed bond and dihedral analysis.

## Requirements
The following Python packages are required to run the script:
- `ase` (Atomic Simulation Environment): For handling crystal structures, I/O operations, and creating supercells.
- `networkx`: For constructing and analyzing graphs based on atomic connectivity.
- `numpy`: For numerical operations.
- `matplotlib`: For visualizing molecular and line graphs.

To install these packages, you can use the command:
```sh
pip install -r requirements.txt
```
or with conda:
```sh
conda env create -f environment.yml
```

## Usage

Put all the `cif` files in `./cifs/` and run:
```sh
python main.py
```
or run:
```sh
python multi_processor.py
```
for multi-processing.
(If there is no `./cifs/`, just create one, or directly run the script, which will create one).



## Notes
- This script filters atoms to include only non-metals based on the provided list of metal elements.
- The supercell transformation matrix is set to `3x3x3` to ensure the proper context for extracting central unit cell molecules.
- Visualization uses `matplotlib` to create both the molecule and line graph plots. Ensure your environment supports graphical output for the plots to display.

