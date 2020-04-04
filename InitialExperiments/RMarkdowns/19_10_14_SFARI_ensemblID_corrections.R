
library(tidyverse)
library(biomaRt)

setwd('./Documents/PhD-Models/FirstPUModel/RMarkdowns')


# Load csv
SFARI_genes = read_csv('./../Data/SFARI/SFARI_genes_08-29-2019.csv')

# Add Ensembl IDs
getinfo = c('hgnc_symbol','ensembl_gene_id','gene_biotype')
mart = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host='feb2014.archive.ensembl.org') ## Gencode v19
datGenes = getBM(attributes=getinfo, filters='hgnc_symbol', values=SFARI_genes$`gene-symbol`, mart=mart)


# Merge SFARI_genes with Ensembl IDs
SFARI_genes = SFARI_genes %>% left_join(datGenes, by=c('gene-symbol'='hgnc_symbol')) %>% mutate('ID' = ensembl_gene_id)


# Convert LRG notation to Ensembl ID
print(paste0(sum(grepl('LRG_', SFARI_genes$ID)), ' genes have LRG references instead of regular Ensembl IDs'))

# Couldn't find a programmatic way to transform LRG to Ensembl IDs ...
# LRG_16 corresponding to the Gene ADA has no ensembl ID, replaced it with the Ensembl ID assigned to the whole gene and not the LRG ID

# 7 LRG IDS found in previous version but not in new version: LRG_176, LRG_289, LRG_334, LRG_340, LRG_472, LRG_512, LRG_7
SFARI_genes$ID = SFARI_genes$ID %>% recode(LRG_111='ENSG00000115415', LRG_130='ENSG00000134982', LRG_138='ENSG00000224389',
                                           LRG_15= 'ENSG00000114062', LRG_16= 'ENSG00000196839', LRG_162='ENSG00000253729',
                                           LRG_176='ENSG00000171316', LRG_196='ENSG00000107099', LRG_199='ENSG00000198947',
                                           LRG_214='ENSG00000196712', LRG_226='ENSG00000184058', LRG_234='ENSG00000104728',
                                           LRG_239='ENSG00000122877', LRG_242='ENSG00000261609', LRG_261='ENSG00000198400',
                                           LRG_289='ENSG00000183873', LRG_293='ENSG00000139618', LRG_307='ENSG00000157404',
                                           LRG_311='ENSG00000171862', LRG_318='ENSG00000141646', LRG_327='ENSG00000145362',
                                           LRG_328='ENSG00000123700', LRG_331='ENSG00000127914', LRG_334='ENSG00000151067',
                                           LRG_340='ENSG00000172893', LRG_343='ENSG00000164362', LRG_355='ENSG00000177663',
                                           LRG_369='ENSG00000169432', LRG_390='ENSG00000198523', LRG_391='ENSG00000155657',
                                           LRG_424='ENSG00000282608', LRG_437='ENSG00000153956', LRG_445='ENSG00000171385',
                                           LRG_472='ENSG00000154767', LRG_493='ENSG00000182197', LRG_512='ENSG00000165671',
                                           LRG_66= 'ENSG00000027697', LRG_7=  'ENSG00000141837', LRG_715='ENSG00000101200',
                                           LRG_8=  'ENSG00000144285', LRG_86='ENSG00000197535')


# 15 Genes didn't return any matches: 
print(paste0('15 Genes didn\'t return any matches: ', paste(sort(SFARI_genes$`gene-symbol`[is.na(SFARI_genes$ID)]), collapse=', ')))

# Couldn't find Ensembl ID for MSNP1AS
SFARI_genes = SFARI_genes %>% mutate(ID=ifelse(`gene-symbol`=='BICDL1', 'ENSG00000135127', 
                                        ifelse(`gene-symbol`=='EMSY',   'ENSG00000158636',
                                        ifelse(`gene-symbol`=='ERBIN',  'ENSG00000112851',
                                        ifelse(`gene-symbol`=='KMT5B',  'ENSG00000110066', 
                                        ifelse(`gene-symbol`=='LNPK',   'ENSG00000144320', 
                                        ifelse(`gene-symbol`=='MIR137', 'ENSG00000284202', 
                                        ifelse(`gene-symbol`=='NEXMIF', 'ENSG00000050030', 
                                        ifelse(`gene-symbol`=='NSMCE3', 'ENSG00000185115',
                                        ifelse(`gene-symbol`=='PATJ',   'ENSG00000132849',
                                        ifelse(`gene-symbol`=='PLPPR4', 'ENSG00000117600',
                                        ifelse(`gene-symbol`=='PRKN',   'ENSG00000185345',
                                        ifelse(`gene-symbol`=='RP11-1407O15.2', 'ENSG00000174093',
                                        ifelse(`gene-symbol`=='TERB2',  'ENSG00000167014',
                                        ifelse(`gene-symbol`=='TSPOAP1','ENSG00000005379', ID)))))))))))))))


write.csv(SFARI_genes, './../Data/SFARI/SFARI_genes_08-29-2019_with_ensembl_IDs.csv', row.names=FALSE)
