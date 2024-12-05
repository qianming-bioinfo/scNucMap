import os
import argparse

parser = argparse.ArgumentParser(description='arguments')
parser.add_argument('--input', '-i', type=str, help='directory of input center file', required=True)
parser.add_argument('--output', '-o', type=str, help='directory of output bed file', required=True)
parser.add_argument('--flank', '-d', type=int, help='the difference of end and start relative to initial center', required=True)
args = parser.parse_args()

try:
    input_dir = args.input
    output_dir = args.output
    flankSize = args.flank
except Exception as e:
    print(e)

# Create output directory if it does not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for file in os.listdir(input_dir):
    if file.endswith(".txt"):
        with open(os.path.join(input_dir, file), "r") as f_in:
            output_file = file.replace(".txt", ".bed")
            with open(os.path.join(output_dir, output_file), "w") as f_out:
                for line in f_in:
                    line = line.strip()
                    chr, bed_center, strand = line.split("\t")
                    newS = str(int(bed_center) - flankSize)
                    newE = str(int(bed_center) + flankSize)
                    f_out.write("\t".join([chr, newS, newE, strand]) + "\n")