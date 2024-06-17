#!/usr/bin/env python3

import argparse
import glob
import os
import random
import subprocess
from natsort import natsorted
import numpy as np
from intervaltree import Interval, IntervalTree

def main():
    parser = argparse.ArgumentParser(description="Generate config file for Indels")
    parser.add_argument('-i', '--indir', type=str, required=True, help='Input BAM directory')
    parser.add_argument('-b', '--bed', type=str, required=True, help='ROI bed')
    parser.add_argument('-g', '--genome', type=str, required=True, help='Reference genome in FASTA format')
    parser.add_argument('--num_variants', type=int, required=True, help='Total number of variants to simulate')
    parser.add_argument('--autosomes', action='store_true', help='Restrict to autosomes')
    args = parser.parse_args()

    if not args.indir or not args.bed or not args.genome or not args.num_variants:
        parser.print_help()
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

    bams = glob.glob(f"{args.indir}/*.bam")
    if not bams:
        print(f"ERROR: no BAM files were found in {args.indir} directory!")
        exit()

    with open(args.bed, 'r') as file:
        lines = file.read().splitlines()

    if args.autosomes:
        lines = [line for line in lines if not (line.startswith('X') or line.startswith('Y'))]

    regions = parse_bed_data(lines)
    variants = generate_variants(regions, args.num_variants, args.genome, samtools)

    for bam in natsorted(bams):
        for variant in variants:
            vaf = np.random.beta(2, 5)
            print(f"{bam}\t{variant[0]}\t{variant[1]}\t{variant[2]}\t{variant[3]}\t{vaf:.4f}")

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
    return ref_seq

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
        return mutated_dup_seq
    else:
        # Generate a random sequence with the specified length
        random_seq = ''.join(random.choice('ACGT') for _ in range(length))
        return random_seq

def generate_indel_size():
    # Generate indel size following a geometric distribution
    p = 0.1  # Probability parameter for the geometric distribution
    size = np.random.geometric(p)
    return min(max(size, 1), 50)  # Ensure size is between 1 and 50

def generate_variants(regions, num_variants, genome, samtools):
    variants = []
    used_intervals = {}
    num_deletions = num_variants // 2
    num_insertions = num_variants - num_deletions

    for region in regions:
        if region[0] not in used_intervals:
            used_intervals[region[0]] = IntervalTree()

    while len(variants) < num_variants:
        region = select_random_region(regions)
        position = random.randint(region[1], region[2])
        indel_size = generate_indel_size()
        interval = Interval(position, position + indel_size)

        if used_intervals[region[0]].overlaps(interval):
            continue

        if len(variants) < num_deletions:
            ref_seq = get_reference_sequence(region[0], position, indel_size, genome, samtools)
            alt_seq = ref_seq[0]  # Only the first base to simulate deletion
            variant = (region[0], position, ref_seq, alt_seq)
        else:
            ref_seq = get_reference_sequence(region[0], position, 1, genome, samtools)
            inserted_seq = generate_inserted_sequence(region[0], position, indel_size, genome, samtools)
            variant = (region[0], position, ref_seq, ref_seq + inserted_seq)

        variants.append(variant)
        used_intervals[region[0]].add(interval)

    return variants

if __name__ == "__main__":
    main()
