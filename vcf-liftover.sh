#!/bin/bash

CHAINFILE=$1
INFILE=$2
OUTFILE=$3
export CHR=$4
STOP_AT_FIRST_CHAIN=$5
METHOD=$6


OURCHAINFILE=$(basename $CHAINFILE)

## Converting chainfile to sensible format

#cat $CHAINFILE | perl -lane 'if($F[0] eq "chain"){if($curchr eq $F[2]){$skip=1}else{$curchr=$F[2];our $current=$F[5];our $coffset=0;$skip=0}} elsif(!$skip) {my $offset=$F[1]-$F[2];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[2];$coffset=$coffset+$offset;}' > $OURCHAINFILE.offset.reverse


if [ "$STOP_AT_FIRST_CHAIN" = "stop" ]; then
	cat $CHAINFILE | perl -lane 'BEGIN{$visited=0;}
	if($F[0] eq "chain"){
		if($visited){exit();}
			if($F[2] ne $F[7] || $F[2] ne $ENV{"CHR"} || $F[4] eq "-" || $F[9] eq "-"){
				$skip=1;
			}
		else{$curchr=$F[2];our $current=$F[5];our $coffset=$F[10]-$F[5];$skip=0;$visited=1;}
		} 
		elsif(!$skip) {my $offset=$F[2]-$F[1];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[1];$coffset=$coffset+$offset;}
		' > $OURCHAINFILE.$CHR.offset
	else
		cat $CHAINFILE | perl -lane '
		if($F[0] eq "chain"){
			if($F[2] ne $F[7] || $F[2] ne $ENV{"CHR"} || $F[4] eq "-" || $F[9] eq "-"){
				$skip=1;
			}
		else{$curchr=$F[2];our $current=$F[5];our $coffset=$F[10]-$F[5];$skip=0;}
		} 
		elsif(!$skip) {my $offset=$F[2]-$F[1];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[1];$coffset=$coffset+$offset;}
		' > $OURCHAINFILE.$CHR.offset
	fi
echo Written $(wc -l $OURCHAINFILE.$CHR.offset) intervals. 
## Parsing the cmdline argument reverse

#if [ "$SENS" = "reverse" ]; then
#	OFFSETS=$OURCHAINFILE.offset.reverse
#	else
OFFSETS=$OURCHAINFILE.$CHR.offset
#fi

#exit()
## Applying the file with Tabix


if [ "$METHOD" =  "sort" ]; then
 echo "Lifting over (method sort)..."
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
 echo "Lifting over (method files)..."
	rm -f liftover.fofnames
	cat $OFFSETS | grep $CHR | while read line;  do 
	(>&2 echo "$line")
	IFS=" " read -r -a fields <<< "$line"
	newstart=$(( ${fields[1]} + ${fields[3]} ))
	newend=$(( ${fields[2]} + ${fields[3]} ))
	closend=$(( ${fields[1]} + 1 ))
	cat <(tabix -h $INFILE fantasy:0-1) <(tabix $INFILE ${fields[0]}:${closend}-${fields[2]} | awk -v offset=${fields[3]} 'BEGIN{OFS="\t"}{$2=$2+offset}1') | bgzip > ${fields[0]}:${newstart}-${newend}.vcf.gz
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

if [ "$METHOD" =  "fast" ]; then
	 echo "Lifting over (method fast)..."
	{
	tabix -H $INFILE
	awk '{print $1, $2, $3, $1, $2+$4, $3+$4, $4}' $OURCHAINFILE.$CHR.offset | sort -k5,5n| awk 'BEGIN{OFS="\t"}{print $1, $5, $6, $2"-"$3":"$7}' | bedtools merge -c 4 -o collapse -i stdin|while read l; do
		j=$(( $j + 1 ))
		a=$(echo $l | tr '\t' ' '| cut -f4 -d' ')
		IFS=',' read -r -a ivals <<< "$a"
		if [ "${#ivals[@]}" -gt 1 ]
			then
			>&2 echo Interval $j sort $a
			for ival in "${ivals[@]}"; do 
				ivl=$(echo $ival | sed 's/:.*//;s/-.*//')
                                ivr=$(echo $ival | sed 's/:.*//;s/.*-//')
				ivl=$(( $ivl + 1 ))
				offset=$(echo $ival | sed 's/.*://')
				tabix $INFILE ${CHR}:${ivl}-$ivr | awk -v offset=$offset 'BEGIN{OFS="\t"}{$3=$2;$2=$2+offset}1'
			done | sort -k2,2n
		else 
			>&2 echo Interval $j nosort $a 
			iv=$(echo ${ivals[0]} | sed 's/:.*//')
			ivl=$(echo $iv | sed 's/:.*//;s/-.*//')
                        ivr=$(echo $iv | sed 's/:.*//;s/.*-//')
                        ivl=$(( $ivl + 1 ))
			offset=$(echo ${ivals[0]} | sed 's/.*://')
			tabix $INFILE ${CHR}:${ivl}-${ivr}  | awk -v offset=$offset 'BEGIN{OFS="\t"}{$3=$2;$2=$2+offset}1' 
		fi 
	done 
	} | bgzip > $OUTFILE
fi
