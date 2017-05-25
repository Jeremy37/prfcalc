#!/bin/bash

DIR=~/src/prfcalc/FANTOM
FANTOM_TPM_FILE=$DIR/hg19.cage_peak_tpm_ann.osc.txt.gz

# The file hg19.cage_peak_tpm_ann.osc.txt.gz.coltpm was prepared by determining the min
# TPM value for each column from the related FANTOM peaks file.
echo -e "Column\tMinTPM\tWeight\tFantomID\tFantomInfo" > hg19.cage_peak_tpm_ann.osc.colweights.txt
cat hg19.cage_peak_tpm_ann.osc.txt.gz.coltpm \
  | awk 'BEGIN {OFS="\t"}{if ($1 ~ /NA/) {print} else {print NR,$1,1/$1,$2}}' \
  | /bin/sed '1,7d' \
  | perl -ane 'chomp; if ($F[3] =~ /.*\.(CNhs[0-9]+)\..*/) { print join("\t", @F[0..2], $1, $F[3])."\n"}' \
   >> hg19.cage_peak_tpm_ann.osc.colweights.txt

# First get the column indexes and weights to use for each epigenome in generating 
# the TSS TPM files.
./getEpigenomeFantomCols.py --roadmap epigenome.metadata.txt \
                            --fantom hg19.cage_peak_tpm_ann.osc.colweights.txt \
                            --output epigenome --debug

# Do the same for enhancers.
./getEpigenomeFantomCols.py --roadmap epigenome.metadata.txt \
                            --fantom hg19.cage_peak_tpm_ann.osc.colweights.txt \
                            --enhancer permissive_enhancers.col.min.tpms.txt \
                            --output epigenome --debug
                            

# This script just runs another python script to get the 
mkdir enh
mkdir tss
submitJobs.py -j getFantomData.all.2.0 --MEM 2000 -o farmOut -c "./getFantomData.allepigenomes.py --tpm 2.0"


# We use three "tiers" of TSS annotations.
# First, we use all TSSes with TPM > 2.0 from FANTOM in a particular cell type
# For any genes with no TSS listed, we next add in the top 3 "true TSSes" determined from
# the FANTOM TSS classifier, irrespective of tissue expression.
# For any genes that still have no TSS listed, we add in all TSSes from Gencode.

cat gencode.v19.tss.txt | awk 'BEGIN {OFS="\t"}{print $1,$2,$3,0}' > gencode.v19.tss.tpm.txt

cd tss
for f in *mintpm.2.0.txt; do
    EID=`echo -n $f | perl -lne '{@a=split(/\./, $_); print $a[0];}'`
    ../mergeFantomTss.pl --fantom $f --ensembl ../truetss.FANTOM.top3.txt > $EID.fantomTSS.mintpm.2.0.plusTop3.txt
    ../mergeFantomTss.pl --fantom $f --ensembl ../gencode.v19.tss.tpm.txt > $EID.fantomTSS.mintpm.2.0.ensg.txt
    hashJoin.pl --scanFile $EID.fantomTSS.mintpm.2.0.ensg.txt --hashFile ../gencode.v19.pc.ids.txt --colHashFile 1 --colScanFile 3 > $EID.fantomTSS.mintpm.2.0.ensg.pc.txt
done
rm *.plusTop3.txt

# Check how many unique genes we have in each file
cat E017.fantomTSS.mintpm.2.0.ensg.pc.txt | cut -f 3 | sort | uniq | wc -l
