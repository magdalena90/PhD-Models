
library(tidyverse)
library(biomaRt)

setwd('./Documents/PhD-Models/FirstPUModel/RMarkdowns')


# Load csv
SFARI_genes = read_csv('./../Data/SFARI/SFARI_genes_01-03-2020.csv')



# Add Ensembl IDs
getinfo = c('hgnc_symbol','ensembl_gene_id','gene_biotype')
mart = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host='feb2014.archive.ensembl.org') ## Gencode v19
datGenes = getBM(attributes=getinfo, filters='hgnc_symbol', values=SFARI_genes$`gene-symbol`, mart=mart)


# Merge SFARI_genes with Ensembl IDs
SFARI_genes = SFARI_genes %>% left_join(datGenes, by=c('gene-symbol'='hgnc_symbol')) %>% mutate('ID'=ensembl_gene_id)


## Convert LRG notation to Ensembl ID: Actually, all the genes with 'LRG_' id have a duplicated entry with their regular ensembl-ID, 
# so there's no need to convert them, we can just remove them
SFARI_genes = SFARI_genes %>% filter(!grepl('LRG_', ID))



# 21 Genes didn't return any matches: 
print(paste0(sum(is.na(SFARI_genes$ID)),' Genes didn\'t return any matches: ', 
             paste(sort(SFARI_genes$`gene-symbol`[is.na(SFARI_genes$ID)]), collapse=', ')))

# So I'm going to copy the Ensembl ID that comes with the dataset, hoping they don't have any alternative IDs, 
# and copy the biomart results that those IDs return
missing_from_biomart = SFARI_genes$`ensembl-id`[is.na(SFARI_genes$ID)]
getinfo = c('ensembl_gene_id','gene_biotype')
mart = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host='feb2014.archive.ensembl.org') ## Gencode v19
datGenes = getBM(attributes=getinfo, filters='ensembl_gene_id', values=missing_from_biomart, mart=mart)

SFARI_genes = SFARI_genes %>% mutate(ID = ifelse(is.na(ID), `ensembl-id`, ID)) %>%
              left_join(datGenes, by=c('ensembl-id'='ensembl_gene_id')) %>%
              mutate(gene_biotype = coalesce(gene_biotype.x, gene_biotype.y)) %>%
              dplyr::select(-gene_biotype.x, -gene_biotype.y)



# Only one gene remaining without gene_biotype and the description says it's a microRNA
SFARI_genes = SFARI_genes %>% mutate(gene_biotype = ifelse(`gene-symbol`=='MIR137', 'microRNA', gene_biotype))



# Write file
write.csv(SFARI_genes, './../Data/SFARI/SFARI_genes_01-03-2020_with_ensembl_IDs.csv', row.names=FALSE)
