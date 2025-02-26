import pysam
import glob
import os
import subprocess
import sys


def sort_bam(input_bam, output_bam=None, n_cpus=1):
    """
    Sort a BAM file.

    :param input_bam: Path to the input BAM file.
    :param output_bam: Path to the output sorted BAM file. If None, replaces the input BAM file name's suffix with '.sorted.bam'.
    :return: Path to the sorted BAM file.
    """
    filedir = os.path.dirname(os.path.abspath(__file__))

    samtools_binary = filedir + "/../samtools/samtools"
    if not os.path.isfile(samtools_binary):
        msg = " ERROR: Missing samtools binary"
        print(msg)
        sys.exit()

    if output_bam is None:
        output_bam = input_bam.rsplit(".", 1)[0] + ".sorted.bam"
    
    # Constructing the samtools sort command
    cmd = [samtools_binary, "sort", "-@", str(n_cpus),  "-T test", "-o", output_bam, input_bam]
    
    try:
        msg = f" INFO: Executing Samtools sorting for {input_bam}"
        print(msg)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        msg = " ERROR: Could not sort bam:", str(e)
        print(msg)

    os.remove(input_bam)
    os.rename(output_bam, input_bam)

    index_bam(input_bam)
    # index_bam(output_bam)


def index_bam(bam_file):
    """
    Indexes a BAM file.

    :param bam_file: Path to the BAM file to be indexed.

    **Example**::

        bam_file = "path/to/your/sorted.bam"
        index_bam(bam_file)
        print("BAM file indexed.")
    """
    pysam.index(bam_file)



def process_all_bam_files(directory, threads):
    """
    Sorts and indexes all BAM files in the specified directory.

    :param directory: The directory to search for BAM files.
    
    **Example**::

        directory = "/path/to/your/directory/with/bamfiles"
        process_all_bam_files(directory)
        print("All BAM files processed.")

    """
    # Fetch all bam files 
    bam_files = glob.glob(os.path.join(directory, "*.bam"))

    for bam_file in bam_files:
        msg = f" INFO: Sorting BAM {bam_file}"
        print(msg)
        sort_bam(input_bam=bam_file, n_cpus=threads)
