import argparse
import os
import sys
import glob
from pybedtools import BedTool

def usage():
    print("Usage: script.py -c <centerBedDir> -r <peakFile> -o <outputDir>")
    print("  -c  Directory containing BED files")
    print("  -r  Path to the sorted open region file in BED format")
    print("  -o  Output directory")
    print("  -h  Show this help message and exit")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Process BED files with pybedtools and intersect with a peak file.")
    parser.add_argument("-c", "--centerBedDir", required=True, help="Directory containing BED files")
    parser.add_argument("-r", "--peakFile", required=True, help="Path to the sorted open region file in BED format")
    parser.add_argument("-o", "--outputDir", required=True, help="Output directory")
    args = parser.parse_args()

    centerBedDir = args.centerBedDir
    peakFile = args.peakFile
    outputDir = args.outputDir

    if not os.path.isdir(centerBedDir):
        print(f"Error: {centerBedDir} is not a directory.")
        usage()
    if not os.path.isfile(peakFile):
        print(f"Error: {peakFile} does not exist or is not a file.")
        usage()

    if not os.path.exists(outputDir):
        os.makedirs(outputDir, exist_ok=True)

    peakPre = os.path.basename(peakFile).split('.')[0]

    for bedFile in glob.glob(os.path.join(centerBedDir, "*.bed")):
        pre = os.path.basename(bedFile).split('.')[0]

        bed = BedTool(bedFile).sort()
        peak = BedTool(peakFile)


        intersected = bed.intersect(peak, wa=True, wb=True, sorted=True)
        
        output_txt = os.path.join(outputDir, f"{pre}_{peakPre}_ovlp.txt")
        with open(output_txt, 'w') as out_f:
            for interval in intersected:
                out_f.write('\t'.join([interval.fields[0], interval.fields[1], interval.fields[3]]) + '\n')

if __name__ == "__main__":
    main()
