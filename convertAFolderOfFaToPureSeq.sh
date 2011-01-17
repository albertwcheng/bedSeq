#!/bin/bash

rm seqsizes.txt

for i in *.fa; do 
convertFaToPureSeq.py $i ${i/.fa/}.seq >> seqsizes.txt
done
