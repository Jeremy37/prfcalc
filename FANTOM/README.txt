The files for each epigenome used in calculating PRF scores are in the subfolders:
tss/
enh/

For TSSes, each file contains, for each gene, all TSSes with FANTOM expression of
2.0 TPM or higher, as well as the TPM level itself. If no FANTOM TSS matched a given
Ensembl gene, or no FANTOM TSS had expression 2 TPM or above, then all Ensembl TSSes are
included in the file with a TPM of 0.0.

For enhancers, each epigenome has a .bed file, where the 4th column is the TPM expression
of the enhancer.

Following is a brief summary of the scripts that produced these files.
get_epigenome_tsses.sh  is the starting point.

The script produces files epigenome.tsscols.txt and epigenome.enhcols.txt, which list the
columns in the main FANTOM files (input/hg19.cage_peak_tpm_ann.osc.txt.gz and
human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz) to use for each
Roadmap epigenome, and the weights to assign to each.
It then extracts these columns for each epigenome to produce files in the tss/ and enh/
folders. Data for TSSes from FANTOM is combined with Ensembl gene annotations to produce
a file with FANTOM gene TPMs followed by Ensembl genes with TPM recorded as 0.0.

prfcalc.py can use one of these tss files as input, and will determine the minimum distance
to any expressed TSS, or alternatively the greater of 5 kb and the minimum distance to any
Ensembl TSS.

