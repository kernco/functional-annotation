#!/bin/bash

rm -rf $1_Jaccard_Temp
mkdir $1_Jaccard_Temp
for bed in Peak_Calls/$1*/*$2; do newbed=${bed/Peak_Calls\//}; ln -f $bed $1_Jaccard_Temp/${newbed//\//_}; done;
cd $1_Jaccard_Temp
parallel --gnu "bedtools jaccard -a {1} -b {2} | awk 'NR>1' | cut -f 3 > {1}.{2}.jaccard" ::: `ls *$1*` ::: `ls *$1*`
ls *$1*.jaccard | xargs grep "" | sed -e s"/\.\///" | perl -pi -e "s/.bed./.bed\t/" | perl -pi -e "s/.jaccard:/\t/" > pairwise.$1.txt
cat pairwise.$1.txt | sed -e "s/$1_//g" | sed -e "s/_macs2_peaks.$2//g" | sed -e 's/\./ /' > pairwise.$1.shortnames.txt
cat pairwise.$1.shortnames.txt | python $(dirname $0)/make-matrix.py > $1.jaccard.matrix
Rscript $(dirname $0)/Jaccard_Graph.R $1.jaccard.matrix ../Results/$1_Peak_Heatmap.png
cd ..
rm -rf $1_Jaccard_Temp
