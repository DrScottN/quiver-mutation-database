import csv
import sys

def main(read_filename, write_filename):
    with open(read_filename, 'r', newline='') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        deduplicatedRecords = dict()

        for row in reader:
            if row["quiver exchange matrix"] not in deduplicatedRecords.keys():
                deduplicatedRecords[tuple(row["quiver exchange matrix"])] = row
            else:
                R = deduplicatedRecords[row["quiver exchange matrix"]]
                for key, value in row.items():
                    if decode(R[key]) == "NULL":
                        R[key] = value
                deduplicatedRecords[row["quiver exchange matrix"]] = R
    
    with open(write_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for quiv, quiv_dict in deduplicatedRecords.items():
            writer.writerow([quiv_dict[key] for key in header])


if __name__ == "__main__":

    args = [j for j  in sys.argv[1:]]
    main(*args)
