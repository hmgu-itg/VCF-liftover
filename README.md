# `vcf-liftover`, a tool to lift over VCFs without going through BED

# TL;DR: 

`./vcf_liftover.sh [chain file] [input.vcf.gz] [output.vcf.gz] [chromosome code]`

* `chain file` : chain file that can be downloaded from the USCS website.
* `input.vcf.gz` : input, tabixed bgzipped VCF format.
* `output.vcf.gz` : output, will be bgzipped VCF.
* `chromosome code` : for example, "chr22". This code has to be identical in both the VCF and the chain file.

## Principle of a `liftOver`

The USCS LiftOver tool allows to translate genomic coordinates from one build to another. For this, it uses a particular file (called a chain file) that describes how positions in both builds correspond to each other.


The [UCSC chain file format](http://genome.ucsc.edu/goldenPath/help/chain.html) can be a bit hard to understand, all the more so that what is written in a chain file is actually a [net](http://genomewiki.ucsc.edu/index.php/Chains_Nets) comprising of several hierarchically arranged chains.

In short, chains are sequences of alignments separated by gaps. So a chain file is easily translatable into a series of genomic intervals with an attached offset. Every SNP located within one such interval will need to be shifted by the corresponding offset. We first write a code snippet that produces a list of intervals and offsets from a chain file.

## Making sense of chain files

The following piece of perl code finds the first chain that describes any given chromosome (as said before, there are many chains per chromosome, and the first, being the one with the highest score, usually spans its entire length). Every ungapped alignment that follows is then translated into an interval.

```bash
cat $CHAINFILE | perl -lane '
if($F[0] eq "chain"){
	if($curchr eq $F[2]){
		$skip=1;
	}else{
		$curchr=$F[2];
		our $current=$F[5];
		our $coffset=$F[10]-$F[5];
		$skip=0;
	}
} elsif(!$skip) {
	my $offset=$F[2]-$F[1];
	print "$curchr $current ", $current+$F[0], " $coffset";
	$current=$current+$F[0]+$F[1];
	$coffset=$coffset+$offset;
}
' > $OURCHAINFILE.offset
```

The code above translates this chain file excerpt:

```
chain 3231099988 chr22 50818468 + 16367188 50806138 chr22 51304566 + 16847850 51244566 23
19744   0       40
36      1       1
```

into the following list of intervals:

```
chr22 16367188 16386932 480662
chr22 16386932 16386968 480702
chr22 16386969 16387000 480702
```

and so forth for all chromosomes. Theoretically, we should descend into further chains instead of using only the first one, as those will further describe mappings within gaps left blank by the topmost chain. 

## Lifting over coordinates in a VCF

What we do next is a bit inefficient: given a chromosome-wide VCF file, we loop over all intervals in the modified chain file, tabixing them out of the VCF and applying the desired offset.

```bash
{ 
tabix -h $INFILE 0:0-1 && cat $OFFSETS | grep $CHR | while read line;  do 
	IFS=" " read -r -a fields <<< "$line"
	tabix $INFILE ${fields[0]}:${fields[1]}-${fields[2]} | awk -v offset=${fields[3]} 'BEGIN{OFS="\t"}{$2=$2+offset;}1'; 
done 
} | bgzip > $OUTFILE
```

## Performance of liftover
### Speed
`vcf-liftover` converts 175k variants in about 10 minutes.

### Accuracy 
UCSC LiftOver also maps positions onto other chromosomes as well. Those typically come from lower-scored chains which `vcf-liftover` ignores. Hence we compare only the cis-mapped positions.

```
~/vcf-liftover/vcf-liftover.sh hg38ToHg19.over.chain 22.vcf.gz 22.liftover.vcf.gz chr22
zgrep -v '^#' 22.vcf.gz | cut -f1,2 | awk '{print $1, $2-1, $2, $1":"$2}' > 22.bed
./liftOver 22.bed hg38ToHg19.over.chain 22.new.bed unmapped.bed
awk '$1=="chr22"' 22.new.bed | cut -f1,3,4 | sed 's/chr22://'> 22.new.formatted.bed
comm -3 <(sort 22.new.2.bed) <(sort 22.vcf.bed) | wc -l
```

We find 593/172026=0.3% of SNP mappings that are unique to `vcf-liftover`, and are presumed to be errors. 34/593=6% of these are unmapped by liftOver. As the reasons for having unmapped SNPs is unknown, it is hard to pinpoint the reason for this. The vast majority of the error SNPs are actually mapped to a different chromosome.