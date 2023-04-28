
# Map SNPs to genes in 35/10Kb window and using NCBI37.3.gene.loc file available from MAGMA website 
magma --annotate window=35,10 \
--snp-loc /path/to/file.snploc \
--gene-loc /path/to/NCBI37.3.gene.loc \
--out /path/to/output/file


# Generate gene-based statistics per chromsome in parallel 
out=/path/to/per_chromsome_output/file
bfile=path/to/plink/bed/bim/fam/files
pval=/path/to/file.pval
gene=path/to/output/file/from/step_1

for i in {1..22}
do
    magma --bfile "$bfile" \
    --pval "$pval" N=1926 \
    --gene-annot "$gene" \
    --batch "$i" chr \
    --gene-model multi=snpwise \
    --out "$out" &
done
wait

# Merge all chr files
magma --merge path/to/folder/containing/per_chromsome/output/files \
--out path/to/final/gene-level/output/file.out
