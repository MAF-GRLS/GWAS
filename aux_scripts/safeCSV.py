import csv
import sys

def process_csv(input_file, output_file):
    with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, quotechar='"', delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer = csv.writer(outfile, delimiter=',', quoting=csv.QUOTE_NONE, escapechar='\\')

        for row in reader:
            new_row = [field.replace(',', ';') for field in row]
            writer.writerow(new_row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_csv(input_file, output_file)

    print("CSV processing complete. Output saved to", output_file)

