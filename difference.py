import csv
import sys

def main(read_filename1, read_filename2, write_filename):
    writeSet = set()
    with open(read_filename1, 'r', newline='') as f:
        reader1 = csv.DictReader(f)
        header1 = reader1.fieldnames
        records1 = set()

        for row in reader1:
            records1.add(row["quiver exchange matrix"])
    with open(read_filename2, 'r', newline='') as f:
        reader2 = csv.DictReader(f)
        header2 = reader2.fieldnames

        for row in reader2:
            if row["quiver exchange matrix"] not in records1:
                writeSet.add(row["quiver exchange matrix"])
    
    with open(write_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for quiv in writeSet:
            writer.writerow(quiv[1:-1].split(' '))


if __name__ == "__main__":

    args = [j for j  in sys.argv[1:]]
    main(*args)