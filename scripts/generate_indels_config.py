#!/usr/bin/env python3

import argparse
import glob
import os
import random
import subprocess
from natsort import natsorted
import numpy as np
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(description="Generate config file for Indels")
    parser.add_argument('-i', '--indir', type=str, help='Input BAM directory')
    parser.add_argument('-l', '--list', type=str, help='Input BAM list file')
    parser.add_argument('-b', '--bed', type=str, required=True, help='ROI bed')
    parser.add_argument('-g', '--genome', type=str, required=True, help='Reference genome in FASTA format')
    parser.add_argument('--num_variants', type=int, required=True, help='Total number of variants to simulate')
    parser.add_argument('--autosomes', action='store_true', help='Restrict to autosomes')
    parser.add_argument('--insertion_rate', type=float, default=0.45, help='Proportion of insertions (default: 0.45)')
    parser.add_argument('--deletion_rate', type=float, default=0.45, help='Proportion of deletions (default: 0.45)')
    parser.add_argument('--delins_rate', type=float, default=0.1, help='Proportion of delins (default: 0.1)')
    parser.add_argument('--chromosome', type=str, help='Specific chromosome to introduce variants')
    parser.add_argument('--vafs', type=str, required=True, help='Comma-separated list of VAFs for different clones')
    parser.add_argument('--std_devs', type=str, required=True, help='Comma-separated list of standard deviations for each clone')
    parser.add_argument('--proportions', type=str, required=True, help='Comma-separated list of proportions for each clone')
    parser.add_argument('--output', type=str, required=True, help='Output PNG file for the VAF distribution plot')
    args = parser.parse_args()

    if not (args.indir or args.list) or not args.bed or not args.genome or not args.num_variants:
        parser.print_help()
        exit()

    total_rate = args.insertion_rate + args.deletion_rate + args.delins_rate
    if total_rate != 1.0:
        print("ERROR: The sum of insertion_rate, deletion_rate, and delins_rate must be 1.0")
        exit()

    if not os.path.exists(args.bed):
        print(f"ERROR: No input ROI bed file found at {args.bed}")
        exit()

    if not os.path.exists(args.genome):
        print(f"ERROR: No reference genome file found at {args.genome}")
        exit()

    samtools = subprocess.getoutput('which samtools')
    if not samtools:
        print("ERROR: samtools was not found on path")
        exit()

    bams = []
    if args.indir:
        bams = glob.glob(f"{args.indir}/*.bam")
        if not bams:
            print(f"ERROR: no BAM files were found in {args.indir} directory!")
            exit()
    elif args.list:
        if not os.path.exists(args.list):
            print(f"ERROR: No input BAM list file found at {args.list}")
            exit()
        with open(args.list, 'r') as file:
            bams = file.read().splitlines()
        if not bams:
            print(f"ERROR: no BAM files were listed in {args.list} file!")
            exit()

    with open(args.bed, 'r') as file:
        lines = file.read().splitlines()

    if args.autosomes:
        lines = [line for line in lines if not (line.startswith('X') or line.startswith('Y'))]

    if args.chromosome:
        lines = [line for line in lines if line.startswith(args.chromosome + "\t")]

    if not lines:
        print(f"ERROR: No regions found in the BED file for the specified chromosome: {args.chromosome}")
        exit()

    regions = parse_bed_data(lines)

    clone_vafs = [float(vaf) for vaf in args.vafs.split(',')]
    std_devs = [float(sd) for sd in args.std_devs.split(',')]
    if len(clone_vafs) != len(std_devs):
        print("ERROR: The number of VAFs must match the number of standard deviations.")
        exit()

    proportions = [float(p) for p in args.proportions.split(',')]
    if len(proportions) != len(clone_vafs):
        print("ERROR: The number of proportions must match the number of clones.")
        exit()
    if not np.isclose(sum(proportions), 1.0):
        print("ERROR: The sum of the proportions must be 1.0")
        exit()

    clone_info = list(zip(clone_vafs, std_devs))
    variants = generate_variants(regions, args.num_variants, args.genome, samtools, args.insertion_rate, args.deletion_rate, args.delins_rate, clone_info, proportions)

    for bam in natsorted(bams):
        for variant in variants:
            print(f"{bam}\t{variant[0]}\t{variant[1]}\t{variant[2]}\t{variant[3]}\t{variant[4]:.4f}")

    plot_vaf_distribution([v[4] for v in variants], args.output, clone_info)

def parse_bed_data(lines):
    regions = []
    for line in lines:
        parts = line.split('\t')
        chr, start, end = parts[0], int(parts[1]), int(parts[2])
        regions.append((chr, start, end))
    return regions

def select_random_region(regions):
    return random.choice(regions)

def get_reference_sequence(chr, position, length, genome, samtools):
    region = f"{chr}:{position}-{position+length-1}"
    cmd = f"{samtools} faidx {genome} {region}"
    ref_seq = subprocess.getoutput(cmd).split('\n')[1]  # Assuming the second line is the sequence
    return ref_seq.upper()

def generate_inserted_sequence(chr, position, length, genome, samtools):
    if random.random() < 1/3:  # 1/3 chance for tandem duplication
        # Select a random nearby region to duplicate
        duplication_start = position - random.randint(1, 1000)
        duplication_start = max(1, duplication_start)  # Ensure the start position is valid
        dup_seq = get_reference_sequence(chr, duplication_start, length, genome, samtools)
        # Introduce random changes with 10% probability
        mutated_dup_seq = ''.join(
            random.choice('ACGT') if random.random() < 0.1 else nt for nt in dup_seq
        )
        return mutated_dup_seq.upper()
    else:
        # Generate a random sequence with the specified length
        random_seq = ''.join(random.choice('ACGT') for _ in range(length))
        return random_seq.upper()

def generate_indel_size():
    # Generate indel size following a geometric distribution
    p = 0.1  # Probability parameter for the geometric distribution
    size = np.random.geometric(p)
    return min(max(size, 1), 50)  # Ensure size is between 1 and 50

def generate_delins(chr, position, genome, samtools):
    while True:
        deletion_size = generate_indel_size()
        if deletion_size > 1:
            break
    insertion_size = random.randint(1, deletion_size - 1)
    ref_seq = get_reference_sequence(chr, position, deletion_size, genome, samtools)
    inserted_seq = generate_inserted_sequence(chr, position, insertion_size, genome, samtools)
    alt_seq = (ref_seq[:1] + inserted_seq).upper()  # Include the first base of the ref_seq followed by the inserted sequence
    return ref_seq, alt_seq

def generate_variants(regions, num_variants, genome, samtools, insertion_rate, deletion_rate, delins_rate, clone_info, proportions):
    variants = []
    used_intervals = {}
    num_insertions = int(num_variants * insertion_rate)
    num_deletions = int(num_variants * deletion_rate)
    num_delins = num_variants - num_insertions - num_deletions

    clone_counts = [int(num_variants * p) for p in proportions]

    for region in regions:
        if region[0] not in used_intervals:
            used_intervals[region[0]] = IntervalTree()

    for idx, clone in enumerate(clone_info):
        clone_vaf, clone_sd = clone
        for _ in range(clone_counts[idx]):
            if len(variants) >= num_variants:
                break
            region = select_random_region(regions)
            position = random.randint(region[1], region[2])
            indel_size = generate_indel_size()
            interval = Interval(position, position + indel_size)

            if used_intervals[region[0]].overlaps(interval):
                continue

            if len(variants) < num_deletions:
                ref_seq = get_reference_sequence(region[0], position, indel_size, genome, samtools)
                alt_seq = ref_seq[0]  # Only the first base to simulate deletion
                variant_type = 'deletion'
            elif len(variants) < num_deletions + num_insertions:
                ref_seq = get_reference_sequence(region[0], position, 1, genome, samtools)
                inserted_seq = generate_inserted_sequence(region[0], position, indel_size, genome, samtools)
                ref_seq = ref_seq.upper()
                alt_seq = (ref_seq + inserted_seq).upper()
                variant_type = 'insertion'
            else:
                ref_seq, alt_seq = generate_delins(region[0], position, genome, samtools)
                variant_type = 'delins'

            vaf = np.random.normal(clone_vaf, clone_sd)
            vaf = min(max(vaf, 0), 1)  # Ensure VAF is within [0, 1]
            variants.append((region[0], position, ref_seq, alt_seq, vaf, variant_type))
            used_intervals[region[0]].add(interval)

    return variants

def plot_vaf_distribution(vafs, output_file, clone_info):
    plt.figure(figsize=(10, 6))
    sns.histplot(vafs, kde=False, bins=40, color='skyblue', label='Simulated VAFs')
    
    for vaf, sd in clone_info:
        sns.kdeplot(np.random.normal(vaf, sd, 20000), label=f'Clone VAF {vaf:.2f}', linewidth=2)
    
    plt.title('VAF Distribution')
    plt.xlabel('Variant Allele Frequency (VAF)')
    plt.ylabel('Count')
    plt.legend()
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    main()
