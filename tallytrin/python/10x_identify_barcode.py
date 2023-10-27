import sys
import logging
import argparse
import cgatcore.iotools as iotools
import pysam

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("10x_identify_barcode.py")

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--whitelist", default=None, type=str, help='Name of whitelist file')
    parser.add_argument("--cmimode", default=0, type=str, help='Run CMI for accuracy evaluation')
    parser.add_argument("--infile", default=None, type=str, help='Nanopore infile FASTQ file')
    parser.add_argument("--outfile", default=None, type=str, help='Name for output FASTQ files')
    parser.add_argument("--barcode", default=16, type=int, help='Length of barcode required')
    parser.add_argument("--umi", default=12, type=int, help='Length of UMI required')
    return parser.parse_args()

def main():
    args = parse_arguments()
    L.info("Arguments:")
    L.info(args)
  
    log = iotools.open_file(f"{args.outfile}.log", "w")
    barcodes = []
    n, y = 0, 0

    with pysam.FastxFile(args.infile) as fh, open(args.outfile, "w") as outfile:
        for record in fh:
            n += 1
            seq = record.sequence
            position = seq.find("AGATCGGAAGAGCGT")
      
            if position != -1:
                umi_start = position - (args.barcode + args.umi)
                umi_end = position - args.barcode
                umi = seq[umi_start:umi_end]
                barcode = seq[position - args.barcode:position]

                if len(umi) == args.umi and len(barcode) == args.barcode:
                    y += 1
                    barcodes.append(barcode)
                    new_seq = seq[:position - (args.barcode + args.umi)]
                    new_quality = record.quality[:position - (args.barcode + args.umi)]
                    outfile.write(f"@{record.name}_{barcode}_{umi}\n{new_seq}\n+\n{new_quality}\n")

    with open(args.whitelist, "w") as out_barcodes:
        for barcode in set(barcodes):
            out_barcodes.write(f"{barcode}\n")

    log.write(f"The number of total reads: {n}\n")
    log.write(f"The number of reads that have a correct barcode and UMI: {y}\n")
    log.write(f"The percentage of total recovered reads is: {((y / n) * 100):.2f}\n")

if __name__ == "__main__":
    main()
