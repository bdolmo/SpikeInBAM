#/home/bdelolmo/Desktop/NA12878/NA12878_S1.bam	chr1	722488	A	ATGC	0.3698

echo '##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##source=RivIndel
##contig=<ID=chrM,length=16571>
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrY,length=59373566>
##INFO=<ID=STATUS,Number=A,Type=String,Description="Germline or Somatic classification">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=FWD,Number=1,Type=Integer,Description="Forward strand supporting reads">
##INFO=<ID=REV,Number=1,Type=Integer,Description="Reverse strand supporting reads">
##INFO=<ID=SB,Number=1,Type=Float,Description="Strand bias (0 to 1)">
##INFO=<ID=KDIV_CONTIG,Number=1,Type=Float,Description="k-mer diversity of the contig">
##INFO=<ID=KDIV_5Flank,Number=1,Type=Float,Description="k-mer diversity of the 5 Upstream sequence">
##INFO=<ID=KDIV_3Flank,Number=1,Type=Float,Description="k-mer diversity of the 3 Downstream sequence">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MBQ,Number=1,Type=Float,Description="Mean Base Quality">
##INFO=<ID=ERR,Number=1,Type=Float,Description="Error rate of the region">
##INFO=<ID=CHIMR,Number=1,Type=Float,Description="Chimeric-read rate of the region">
##INFO=<ID=SCR,Number=1,Type=Float,Description="Soft-clipped rate of the region">
##INFO=<ID=GC,Number=1,Type=Float,Description="GC-content of the contig">
##INFO=<ID=HOML,Number=1,Type=Integer,Description="Longest homopolymer run">
##INFO=<ID=HOMN,Number=1,Type=Integer,Description="Number of homopolymers">
##INFO=<ID=DNTDL,Number=1,Type=Integer,Description="Longest di-nucleotide run">
##INFO=<ID=DNTDN,Number=1,Type=Integer,Description="Number of di-nucleotide runs">
##INFO=<ID=SOURCE,Number=1,Type=String,Description="Caller name">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of Coverage of the variant">
##bcftools_normVersion=1.13+htslib-1.19.1-19-g5d2c3f72
##bcftools_normCommand=norm -f /home/bdelolmo/REF_DIR/hg19/ucsc.hg19.fasta -o test.somatic.rivindel.vcf -O v test.somatic.rivindel.vcf.final.vcf; Date=Mon Aug 12 14:28:02 2024
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878' > test.vcf

cut -f 2,3,4,5,6 somatic.indels.config | awk '{ print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tAF="$5"\t.\t."}' >> test.vcf

bcftools norm -f ~/REF_DIR/hg19/ucsc.hg19.fasta -o known.somatic.norm.vcf -O v  test.vcf -cs

