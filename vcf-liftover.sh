#!/bin/bash

CHAINFILE=$1
INFILE=$2
OUTFILE=$3
CHR=$4
#SENS=$5

OURCHAINFILE=$(basename $CHAINFILE)

## Converting chainfile to sensible format

#cat $CHAINFILE | perl -lane 'if($F[0] eq "chain"){if($curchr eq $F[2]){$skip=1}else{$curchr=$F[2];our $current=$F[5];our $coffset=0;$skip=0}} elsif(!$skip) {my $offset=$F[1]-$F[2];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[2];$coffset=$coffset+$offset;}' > $OURCHAINFILE.offset.reverse
echo H1
cat $CHAINFILE | perl -lane 'if($F[0] eq "chain"){if($curchr eq $F[2]){$skip=1}else{$curchr=$F[2];our $current=$F[5];our $coffset=$F[10]-$F[5];$skip=0}} elsif(!$skip) {my $offset=$F[2]-$F[1];print "$curchr $current ", $current+$F[0], " $coffset"; $current=$current+$F[0]+$F[1];$coffset=$coffset+$offset;}' > $OURCHAINFILE.offset
echo H2
## Parsing the cmdline argument reverse

#if [ "$SENS" = "reverse" ]; then
#	OFFSETS=$OURCHAINFILE.offset.reverse
#	else
	OFFSETS=$OURCHAINFILE.offset
#fi


## Applying the file with Tabix

{ 
tabix -h $INFILE 0:0-1 && cat $OFFSETS | grep $CHR | while read line;  do 
	(>&2 echo "$line")
	IFS=" " read -r -a fields <<< "$line"
	tabix $INFILE ${fields[0]}:${fields[1]}-${fields[2]} | awk -v offset=${fields[3]} 'BEGIN{OFS="\t"}{$3=$2;$2=$2+offset}1'; 
done 
} | bgzip > $OUTFILE