import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

input_file = "./data/filename-smi.csv"
output_file = "./data/uni-smi.csv"


def read_smiles(file):
    smiles_to_filenames = defaultdict(list)  # Dictionary to store SMILES as keys and filenames as values
    with open(file, mode="r") as infile:
        reader = csv.reader(infile)
        next(reader)  # Skip the header row
        for row in reader:
            filename, smiles = row
            smiles_to_filenames[smiles].append(filename)
    return smiles_to_filenames


def write_smiles(file, smiles_map):
    with open(file, mode="w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["SMILES", "Filenames"])  # Write the header row
        for smiles, filenames in smiles_map.items():
            writer.writerow([smiles, ";".join(filenames)])  # Join filenames with a semicolon


def main():
    with ThreadPoolExecutor(max_workers=16) as executor:
        future_read = executor.submit(read_smiles, input_file)
        smiles_to_filenames = future_read.result()
        executor.submit(write_smiles, output_file, smiles_to_filenames)


if __name__ == "__main__":
    main()
