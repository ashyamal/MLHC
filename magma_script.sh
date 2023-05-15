#Script for running magma gene annotation, gene analysis, and gene-set analysis. 
#Input textfiles were formatted either via the Python code in the repo or awk Unix commands (not shown here for brevity)

magma --annotate --snp-loc SNP_Loc.txt --gene-loc g1000_eur --out SNP_to_gene
magma --bfile g1000_eur --pval SNP_p_vals.txt N=3866 --gene-annot SNP_to_gene.genes.annot --out brain_magma
magma --gene-results brain_magma.genes.raw --set-annot Tx_subtypes_grouped.txt --out image_tx_subtypes
