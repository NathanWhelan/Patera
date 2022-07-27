#!/bin/bash

###Commands to build phylogenetic datasets from STACKS###

#This script assumes you did standard filtering with populations, then made a popmap where all individuals had a seperate population
#designation and made a whitelist from the standard filtering. This is required because STACKS only outputs phylips by populaion.
#There may be a better way to do this that I don't know about.

##Code to create whitelist: grep -v "#" populations.snps.vcf |awk -F "\t" '{print $1}' |sort |uniq >whitelist.txt

##Then run populations with new popmap, -r 1.0, and whitelist.

###RUN THIS SCRIPT FROM INSIDE THE POPULATIONS OUTPUT FOLDER

##Must be in PATH: AMAS.py, singleline.pl, noGaps_nucleotides.sh, iqtree, GNU parallel

#First remove last line of phylip
mkdir phylogenetics
mv *.phylip phylogenetics
cd phylogenetics
conda activate
AMAS.py convert -i populations.all.phylip -f phylip-int -d dna -c 10 -u phylip
sed 's/DNA, //' populations.all.partitions.phylip >populations.all.partitions.simple.txt
AMAS.py split -l populations.all.partitions.simple.txt -u fasta -i populations.all.phylip-out.phy -f phylip -d dna
for FILE in *.fas; do singleline.pl $FILE >$FILE.sl; done
rename .fas.sl .fas *.fas.sl
noGaps_nucleotides.sh
ls *.fas | parallel -j 37 'iqtree -s {} -bb 1000 -nt 2 -m MFP --runs 7 -pers 0.2 -nstop 500'
