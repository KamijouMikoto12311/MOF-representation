# README

## Overview
This script (`main.py`) processes crystal structure data, generates a supercell, identifies atoms within a central unit cell, and creates molecular subgraphs based on atomic connectivity. It extracts the ligands from the crystal structure and visualizes both the molecule and line graphs with detailed bond analysis.

## Requirements
The following Python packages are required to run the script:
- `ase` (Atomic Simulation Environment): For handling crystal structures, I/O operations, and creating supercells.
- `networkx`: For constructing and analyzing graphs based on atomic connectivity.
- `numpy`: For numerical operations.
- `matplotlib`: For visualizing molecular and line graphs.

To install these packages, you can use the command:
```sh
pip install ase networkx numpy matplotlib
```

## Usage
### Input
The script reads the input crystal structure from a `.cif` file named `3D-solvent.cif` by defaault. Change this into the cif file you want to analyze. Ensure that the input file is located in the same directory as the script or modify the code accordingly.

### Running the Script
To run the script:
```sh
python main.py
```
Ensure the necessary input file (`3D-solvent.cif`) is in the correct location.

### Output
- Extracted molecules are saved in the `ligands_xyz` directory with filenames like `extracted_molecule_X.xyz`, where `X` is the index of the molecule.
- Visualization of the molecules and their line graphs is displayed using `matplotlib`.


## Notes
- This script filters atoms to include only non-metals based on the provided list of metal elements.
- The supercell transformation matrix is set to `3x3x3` to ensure the proper context for extracting central unit cell molecules.
- Visualization uses `matplotlib` to create both the molecule and line graph plots. Ensure your environment supports graphical output for the plots to display.

