import sys
import os
import csv

def read_lhloglikelihood(file_path):
    lls = []
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            lls.append(float(row["LHLogLikelihood"]))
    return lls


def read_a_specific_line_as_float(filename):
    try:
        with open(filename, 'r') as f:
            for i, line in enumerate(f, start=1):
                if i == 33:
                    return float(line.strip())
    except Exception as e:
        print(f"Error reading {filename}: {e}")
    return None


def read_a_folder(folder, start_index, end_index):
    lws = []
    for i in range(start_index, end_index + 1):
        filepath = os.path.join(folder, f"{i}.out")
        value = read_a_specific_line_as_float(filepath)
        if value is not None:
            lws.append(value)
        else:
            lws.append(float('-inf'))
    return lws


def sort_indices(values):
    return sorted(range(len(values)), key=lambda i: values[i], reverse=True)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python process_lls.py <log> <folder> <start_index> <end_index>")
        sys.exit(1)

    file_path = sys.argv[1]
    ism = read_lhloglikelihood(file_path)
    print("Order under the ISM: " + str(sort_indices(ism)))
    lws = read_a_folder(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    dsm = [ism[i] + lws[i] for i in range(len(ism))]
    print("Order under the DSM: " + str(sort_indices(dsm)))
