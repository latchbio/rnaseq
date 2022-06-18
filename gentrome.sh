#! /bin/bash

genome_ref=$1
transcipt_ref=$2

grep "^>" $genome_ref | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat $transcript_ref $genome_ref  > /root/gentrome.fa