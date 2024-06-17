# SpikeInBAM - Fast Variant Simulation on BAMs by Spike-In

SpikeInBAM is a fast tool to simulate variants on already available BAMs.

## Requirements and installation:

- You will need a genome reference in FASTA format.
- Install the required python libraries:
```
    pip3 install -r requirements.txt
```

## Usage:

You will need a config file referred as `variants.txt` specifying the variants to be simulated.
There should be a variant per line. Here an example:
```
    /home/user/input_folder/sample1.bam  chr7   55242467   GAATTAAGAGAAGCAACA   GTTGCT
    /home/user/input_folder/sample2.bam  chr18   29104352   A   T
    /home/user/input_folder/sample3.bam  chr1   156082722    156089883    DEL    SOS1_7_to_8    1   multiple
    /home/user/input_folder/sample4.bam  chr7   150652435    150661275    DUP	TAZ_1_to_4	1	multiple
```

To execute the program:

```
    python3 spikeinbam.py  --variants variants.txt --reference /path/to/ref.fasta --threads 4 --suffix .simulated --output <output_directory>
```

The suffix parameter by default is defined as ".simulated" just before ".bam" (e.g sample1.simulated.bam).

All bam files will be placed at outdir and will display the suffix *.simulated.bam


## Additional scripts:

You may consider useful to generated config files scripts/generate_config.py.
It will create a config file supporting randomly placed single-exon or multiple-exon CNVs (--mode param).

```
    python3 generate_cnv_config.py --indir /path/to/directory_with_bam_files -c <del/dup> -b  /path/to/gene_panel.bed  --mode <single/multiple> > output.config
```

## Motivation:

There are already several programs that can generate synthetic datasets.
However, not all of them can be used to "inject" variants over existing BAMs. This feature is desirable when testing the specificity of any bioinformatic tool.
In addition, some of these tools tend to include some sophisticated steps (assembly, re-mapping, etc) that are computationally expensive.

SpikeInBAM is an easy to use tool to simulate germline variants in a fast manner.

SpikeInBAM can:
- Simulate SNV
- Simulate Indel and complex Indels
- Simulate CNVs. It modifies the VAFs of any SNV that overlaps the CNV to be simulated

## Limitations:

- Simulate Structural Variants (SVs) supported by sequenced breakpoints.
- Somatic variants.
- Copy-number changes greater than 3. No homo/hemizygous deletions are supported yet.
