import os
import csv
import subprocess
import shutil
from concurrent.futures import ThreadPoolExecutor


# extract and rename xyz files
# re-annotate xyz files to the file name
# generate csv (filename, smi)


def process_file(file, root, base_dir, output_dir):
    if file.endswith(".xyz"):
        relative_path = os.path.relpath(root, base_dir)
        dir_parts = relative_path.split(os.sep)
        if len(dir_parts) >= 2:
            new_name = f"{dir_parts[-2]}-{dir_parts[-1]}.xyz"
        else:
            new_name = f"{dir_parts[0]}-{file}"

        src = os.path.join(root, file)
        dst = os.path.join(output_dir, new_name)
        shutil.copy(src, dst)
        print(f"Copied and renamed: {src} -> {dst}")


def extract_and_rename_xyz_files(base_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = []
        for root, _, files in os.walk(base_dir):
            for file in files:
                futures.append(executor.submit(process_file, file, root, base_dir, output_dir))
        for future in futures:
            future.result()


def re_annotate_xyz_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".xyz"):
            filepath = os.path.join(directory, filename)
            with open(filepath, "r") as file:
                lines = file.readlines()

            lines[1] = filename + "\n"

            with open(filepath, "w") as file:
                file.writelines(lines)


if __name__ == "__main__":
    base_dir = "./data/ligands"
    output_dir = "./data/ligands-xyz"
    output_csv = "./data/filename-smi.csv"

    extract_and_rename_xyz_files(base_dir, output_dir)
    re_annotate_xyz_files(output_dir)

    with open(output_csv, mode="w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Filename", "SMILES"])

        for filename in os.listdir(output_dir):
            if filename.endswith(".xyz"):
                input_file = os.path.join(output_dir, filename)

                try:
                    command = ["obabel", input_file, "-osmi", "--canonical"]
                    result = subprocess.run(command, capture_output=True, text=True)

                    if result.returncode != 0:
                        print(f"Error converting {filename}: {result.stderr.strip()}")
                        continue

                    smiles = result.stdout.strip().split("\t")[0]

                    csv_writer.writerow([filename, smiles])
                    print(f"Successfully converted {filename} to SMILES and saved to CSV.")

                except Exception as e:
                    print(f"An unexpected error occurred with {filename}: {str(e)}")
