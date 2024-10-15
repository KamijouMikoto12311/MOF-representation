# README

## Overview
This script (`super.py`) processes crystal structure data, generates a supercell, identifies atoms within a central unit cell, and creates molecular subgraphs based on atomic connectivity. It extracts the ligands from the crystal structure and visualizes both the molecule and line graphs with detailed bond analysis.

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
The script reads the input crystal structure from a `.cif` file named `str_m2_o2_o11_pcu_sym.25.cif` by defaault. Change this into the cif file you want to analyse. Ensure that the input file is located in the same directory as the script or modify the code accordingly.

### Steps Performed
1. **Load Structure**: Reads the crystal structure using ASE.
2. **Non-Metal Extraction**: Filters out non-metal atoms from the structure.
3. **Supercell Creation**: Generates a `3x3x3` supercell of the crystal.
4. **Central Unit Cell Identification**: Identifies atoms in the central unit cell of the supercell.
5. **Graph Construction**: Builds a graph where nodes are atoms and edges represent atomic bonds.
6. **Subgraph Extraction**: Extracts connected components (molecules) from the central unit cell.
7. **Line Graph Analysis**: Generates a line graph from the molecular subgraph to analyze bond angles.
8. **Visualization**: Visualizes extracted molecules and their corresponding line graphs.
9. **Output**: Extracted molecules are saved in `.xyz` format in the `ligands_xyz` directory.

### Running the Script
To run the script:
```sh
python main.py
```
Ensure the necessary input file (`str_m2_o2_o11_pcu_sym.25.cif`) is in the correct location.

### Output
- Extracted molecules are saved in the `ligands_xyz` directory with filenames like `extracted_molecule_X.xyz`, where `X` is the index of the molecule.
- Visualization of the molecules and their line graphs is displayed using `matplotlib`.

## Utility Functions
- `utils/visualize_molecule_graph.py`: Contains functions for visualizing molecular graphs.
- `utils/LineGraph.py`: Contains functions for creating and visualizing line graphs of subgraphs, including bond angles.
- `utils/compare.py`: Includes functionality for removing duplicate subgraphs.

## Notes
- This script filters atoms to include only non-metals based on the provided list of metal elements.
- The supercell transformation matrix is set to `3x3x3` to ensure the proper context for extracting central unit cell molecules.
- Visualization uses `matplotlib` to create both the molecule and line graph plots. Ensure your environment supports graphical output for the plots to display.

