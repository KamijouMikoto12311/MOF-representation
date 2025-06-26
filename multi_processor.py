import os
import warnings
from functools import partial  # For future use
from concurrent.futures import ProcessPoolExecutor
from utils.process_cif import process_cif

warnings.filterwarnings("ignore")


def extract_ligands(input_dir, linegraph=False):
    file_paths = [
        os.path.join(input_dir, filename)
        for filename in os.listdir(input_dir)
        if os.path.isfile(os.path.join(input_dir, filename))
    ]

    process_cif_L = partial(process_cif, linegraph=linegraph)

    with ProcessPoolExecutor(max_workers=32) as executor:
        executor.map(process_cif_L, file_paths, chunksize=1)


if __name__ == "__main__":
    input_dir = "./data/cifs"
    extract_ligands(input_dir=input_dir, linegraph=True)
