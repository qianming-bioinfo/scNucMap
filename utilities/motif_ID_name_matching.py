import argparse
from pyjaspar import jaspardb

def parse_args():
    parser = argparse.ArgumentParser(description="Query motif names from JASPAR database using motif IDs.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input file containing motif IDs.")
    parser.add_argument('-o', '--output', required=True, help="Path to the output file to save the results.")
    return parser.parse_args()

def main():
    args = parse_args()
    
    input_file = args.input
    output_file = args.output

    jdb_obj = jaspardb(release='JASPAR2022')

    with open(output_file, 'w') as output:
        
        with open(input_file, 'r') as input:
            for line in input:
                motif_id = line.strip()
                try:
                    motif = jdb_obj.fetch_motif_by_id(motif_id)
                    motif_name = motif.name
                    output.write(f"{motif_id}\t{motif_name}\n")
                except Exception as e:
                    print(f"Error fetching motif for ID {motif_id}: {e}")
                    output.write(f"{motif_id}\tNot Found\n")

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()
