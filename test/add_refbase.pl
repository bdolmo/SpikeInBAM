#!/usr/bin/env perl
use strict;

my $infile = $ARGV[0];

open IN, "<", $infile;
while (my $line=<IN>) {
    #ENST00000275493:E746_P753delinsLS	7	55242466	GAATTAAGAGAAGCAACATCTC	CTCT
    my @tmp = split("\t", $line);

    my $pos = $tmp[2]-1;
    my $coordinates = "$tmp[1]:$pos-$pos";

    my $base = `samtools faidx $coordinates`;
    chomp $base;

    print "$tmp[0]\t$tmp[1]\t$pos\t$base$tmp[3]\t$base$tmp[4]\n";

}
close IN;