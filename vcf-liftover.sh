#!/bin/bash

CHAINFILE=$1
INFILE=$2
OUTFILE=$3
export CHR=$4
STOP_AT_FIRST_CHAIN=$5
METHOD=$6

#PICARD="gsub 10g -I /software/jre1.8.0_74/bin/java -Xmx10000m -jar /lustre/scratch115/projects/t144_helic_15x/analysis/HA/vs_1x/picard.jar"
PICARD="/software/jdk1.8.0_74/bin/java -Xmx30g -jar /lustre/scratch115/projects/t144_helic_15x/analysis/HA/vs_1x/picard.jar"
OURCHAINFILE=$(basename $CHAINFILE)

## Converting chainfile to sensible format

#cat $CHAINFILE | perl -lane 'if($F[0] eq "chain"){if($curchr eq $F[2]){$skip=1}else{$curchr=$F[2];our $current=$F[5];our $coffset=0;$skip=0}} elsif(!$skip) {my $offset=$F[1]-$F[2];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[2];$coffset=$coffset+$offset;}' > $OURCHAINFILE.offset.reverse


if [ "$STOP_AT_FIRST_CHAIN" = "stop" ]; then
	cat $CHAINFILE | perl -lane 'BEGIN{$visited=0;}
	if($F[0] eq "chain"){
		if($visited){exit();}
			if($F[2] ne $F[7] || $F[2] ne $ENV{"CHR"}){
				$skip=1;
			}
		else{$curchr=$F[2];our $current=$F[5];our $coffset=$F[10]-$F[5];$skip=0;$visited=1;}
		} 
		elsif(!$skip) {my $offset=$F[2]-$F[1];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[1];$coffset=$coffset+$offset;}
		' > $OURCHAINFILE.offset
	else
		cat $CHAINFILE | perl -lane '
		if($F[0] eq "chain"){
			if($F[2] ne $F[7] || $F[2] ne $ENV{"CHR"}){
				$skip=1;
			}
		else{$curchr=$F[2];our $current=$F[5];our $coffset=$F[10]-$F[5];$skip=0;}
		} 
		elsif(!$skip) {my $offset=$F[2]-$F[1];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[1];$coffset=$coffset+$offset;}
		' > $OURCHAINFILE.offset
	fi
echo Written $(wc -l $OURCHAINFILE.offset) intervals. Lifting over...
## Parsing the cmdline argument reverse

#if [ "$SENS" = "reverse" ]; then
#	OFFSETS=$OURCHAINFILE.offset.reverse
#	else
OFFSETS=$OURCHAINFILE.offset
#fi

#exit()
## Applying the file with Tabix


if [ "$METHOD" =  "sort" ]; then
 { 
tabix -h $INFILE 0:0-1 && cat $OFFSETS | grep $CHR | while read line;  do 
	(>&2 echo "$line")
	IFS=" " read -r -a fields <<< "$line"
	closend=$(( ${fields[1]} + 1 ))
	tabix $INFILE ${fields[0]}:${closend}-${fields[2]} | awk -v offset=${fields[3]} 'BEGIN{OFS="\t"}{$3=$2;$2=$2+offset}1'; 
done 
} | bgzip > $OUTFILE
echo "Done. Sorting..."
$PICARD SortVcf I=$OUTFILE O=$OUTFILE.sorted.vcf.gz CREATE_INDEX=false
mv $OUTFILE.sorted.vcf.gz $OUTFILE
fi


if [ "$METHOD" =  "files" ]; then
	rm -f liftover.fofnames
	cat $OFFSETS | grep $CHR | while read line;  do 
	(>&2 echo "$line")
	IFS=" " read -r -a fields <<< "$line"
	newstart=$(( ${fields[1]} + ${fields[3]} ))
	newend=$(( ${fields[2]} + ${fields[3]} ))
	closend=$(( ${fields[1]} + 1 ))
	cat <(tabix -h $INFILE fantasy:0-1) <(tabix $INFILE ${fields[0]}:${closend}-${fields[2]} | awk -v offset=${fields[3]} 'BEGIN{OFS="\t"}{$3=$2;$2=$2+offset}1') | bgzip > ${fields[0]}:${newstart}-${newend}.vcf.gz
	nsnps=$(zgrep -c -v '^#' ${fields[0]}:${newstart}-${newend}.vcf.gz)
	if [ $nsnps -lt 1 ]; then
		echo "No SNPs in ${fields[0]}:${newstart}-${newend}"
		rm ${fields[0]}:${newstart}-${newend}.vcf.gz

	else
		tabix -f -p vcf ${fields[0]}:${newstart}-${newend}.vcf.gz
		echo $newstart ${fields[0]}:${newstart}-${newend}.vcf.gz >> liftover.fofnames
	fi
done

echo "Done. Sorting..."
sort -k1,1n liftover.fofnames | cut -d' '  -f2 | sponge liftover.fofnames
bcftools concat -f liftover.fofnames -a -Oz -o surprise.vcf.gz
cat liftover.fofnames | while read f; do rm $f; done
rm liftover.fofnames
rm $CHR:*.tbi
fi

