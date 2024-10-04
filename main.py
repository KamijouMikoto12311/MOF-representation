from ase.io import read, write
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms
import numpy as np

# Step 1: Read the MOF CIF file
structure = read("str_m2_o2_o11_pcu_sym.25.cif")

# Step 2: Filter out metal atoms, leaving only organic atoms (C, H, N, O)
organic_atoms = [
    atom.index
    for atom in structure
    if atom.symbol not in ["Zn", "Cu", "Fe", "Al", "Ti", "Zr", "Co", "Ni"]
]

# Create a new Atoms object containing only the organic atoms
organic_structure = structure[organic_atoms]

# Step 3: Generate cutoff radii for the organic atoms
ori_cutoffs = natural_cutoffs(organic_structure)
cutoffs = [x + 0.4 for x in ori_cutoffs]

# Step 4: Create a NeighborList considering PBC
neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
neighbor_list.update(organic_structure)

# Step 5: Find all connected clusters of organic atoms (linkers)
linkers = []
visited = set()

for atom_idx in range(len(organic_structure)):
    if atom_idx not in visited:
        # Find neighbors of the current atom, including across periodic boundaries
        indices, offsets = neighbor_list.get_neighbors(atom_idx)
        linker = [(atom_idx, np.zeros(3))]  # Store index and offset (to apply later)
        visited.add(atom_idx)

        # Recursively add neighbors to the linker group
        for neighbor_idx, offset in zip(indices, offsets):
            if neighbor_idx not in visited:
                linker.append((neighbor_idx, offset))
                visited.add(neighbor_idx)

        linkers.append(linker)

# Step 6: Output each organic linker as a separate structure
for i, linker_data in enumerate(linkers):
    # Extract indices and offsets for the linker
    linker_indices, offsets = zip(*linker_data)

    # Extract the atoms corresponding to the linker
    linker_atoms = organic_structure[list(linker_indices)].copy()

    # Step 7: Adjust atom positions to move them together into the same unit cell
    for j, offset in enumerate(offsets):
        linker_atoms.positions[j] += np.dot(
            offset, organic_structure.get_cell()
        )  # Apply periodic offset

    # Ensure that the adjusted structure uses PBC
    linker_atoms.set_pbc(True)

    # Write the adjusted linker to an output file
    output_filename = f"organic_linker_{i+1}.xyz"
    write(output_filename, linker_atoms)
    print(f"Outputted organic linker {i+1} to {output_filename}")
