import argparse
import subprocess
import os
import time
from modules.utils import process_all_bam_files

def main():
    parser = argparse.ArgumentParser(description='SpikeInBAM: simulate variants (spike in) to already available BAMs')
    parser.add_argument('--variants', '-v', required=True, help='Variants to be simulated')
    parser.add_argument('--threads', '-t', default=1, help='Total number of threads')
    parser.add_argument('--reference', '-r', required=True, help='Genome reference in FASTA format')
    parser.add_argument('--output', '-o', required=True, help='Output directory where all simulated BAMs will be placed')
    parser.add_argument('--suffix', '-s', required=True, help='Suffix')

    args = parser.parse_args()
    start_time = time.time()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    spikeinbam_exe = os.path.join(script_dir, "src", "spikeinbam")

    command = [
        spikeinbam_exe,
        args.variants,
        args.reference,
        args.output,
        args.suffix,
    ]
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    
    try:
        msg = " INFO: Executing SpikeInBAM:", " ".join(command)
        print(msg)
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        msg = " ERROR: Could not execute SpikeInBAM:", str(e)
        print(msg)

    process_all_bam_files(args.output, args.threads)

    end_time = time.time()  # Record the end time
    total_time = end_time - start_time  # Calculate the total execution time
    msg = f" INFO: Total execution time: {total_time:.2f} seconds"
    print(msg)


if __name__ == "__main__":
    main()
