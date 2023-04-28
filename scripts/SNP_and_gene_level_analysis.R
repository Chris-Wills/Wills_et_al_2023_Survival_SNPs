library(data.table)
library(tidyverse)
library(gwasurvivr)
library(survSNP)

setwd("path/to/working/directory/")

#### GWAS for SNP-level analysis ####
# Read COIN/COIN-B patient data
data <- fread("path/to/patient/data.tsv")

# Perform multivariate GWAS of all SNPs using gwasurvivr package
# include 11 prognostic covariates (tumour location encoded as 7 binary dummy variables)
# minor allele frequency filter at 0.01
# 1926 samples included after filtering for complete covariate data records
options("gwasurvivr.cores"=16)
plinkCoxSurv(bed.file = "path/to/plink/bed/bim/fam/files",
             covariate.file = data,
             id.column = "ID",
             event = "EVENT",
             time.to.event = "SURVIVAL", 
             covariates = c("PS", "RS", "WBC", "PLT", "ALKP", "NMS", "MLIV", "RC",
                            "RJ", "SC", "R", "MG", "LC", "TC", "SA", "METSTIME", "METSCAT"),
             maf.filter = 0.01,
             out.file = "path/to/out/file")


# Read results
gwas <- fread("path/to/out/file")

# Read list of 181 SNPs from Fernandez-Rozadilla et al. 2023 that passed QC in COIN/COIN-B
snps <- fread("path/to/risk/snps/list.tsv")

# filter gwas results for selected snps and write results
gwas %>%
  subset(RSID %in% snps$SNP) %>%
  fwrite("path/to/final/results/file",
         col.names = T,
         row.names = F,
         quote = F,
         sep = "\t")


#### Prepare summary statistics file for MAGMA gene-level analysis ####
# create snploc file for magma 
gwas[,c("RSID", "CHR", "POS")] %>% 
  fwrite( "path/to/file.snploc",
          col.names = F,
          row.names = F, sep = "\t",
          quote = F)

#create pval file for magma 
gwas[,c("RSID", "PVALUE"),] %>% 
  fwrite("path/to/file.pval",
       col.names = F,
       row.names = F, sep = "\t", quote = F)

### RUN MAGMA shell file 

# load .out file
all_gene <- fread("path/to/MAGMA/output.genes.out")
gene_loc <- fread("path/to/NCBI37.3.gene.loc")

# Add gene symbols using the NCBI37.3.gene.loc file and write to file
all_gene <- mutate(all_gene, GENE_NAME = gene_loc$V6[match(all_gene$GENE, gene_loc$V1)])
fwrite(all_gene,
       "path/to/final/gene/results/output/file.tsv",
       col.names = T,
       row.names = F,
       quote = F,
       sep = "\t")

# Subset genes from Fernandez-Rozadilla et al. 2023 (52 out of 53 available) and write to file
genes <- fread("path/to/risk/genes")
all_gene %>%
  subset(all_gene, GENE_NAME %in% genes$gene) %>%
  fwrite("path/to/final/risk/gene/results/output/file.tsv",
         col.names = T,
         row.names = F,
         quote = F,
         sep = "\t") 




#### Power calculations for rs10161980 ####
sim.snp.expsurv.power(GHR = 1.24,
                      n = 1926,
                      raf = 0.374,
                      erate = 1435/1926,
                      B = 0,
                      pilm=0.5,
                      lm=1,
                      model="recessive",
                      test="recessive",
                      alpha=0.05/233)





