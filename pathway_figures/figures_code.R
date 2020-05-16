setwd("/media/barbitoff/DATA/Working issues/WES/Pregnancy_GWAS/huge")

library(ggplot2)
library(reshape2)
library(ggforce)
library(VennDiagram)
library(colorRamps)
library(RColorBrewer)

jet = matlab.like2(100)
mypal <- colorRampPalette(c('white', 'red'))(100)
mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)

# Creation of the tables
# for i in MSIGDB_HUGE_* ; do grep -P '^Gene\ Set|^[A-Z][A-Z]+' $i | sed 's/ /_/g' | sed 's/\#/num/g' > huge/${i%%.tsv}.reformatted.tsv ; done

# Visualizing MSigDB online enrichment data
gest = read.table('MSIGDB_HUGE_Gestational_Diabetes.reformatted.tsv', sep='\t',
                  header=T, quote="")[1:5, ]

preterm = read.table('MSIGDB_HUGE_Preterm_Birth.reformatted.tsv', sep='\t',
                  header=T, quote="")[1:5, ]

eclam = read.table('MSIGDB_HUGE_Preeclampsia.reformatted.tsv', sep='\t',
                  header=T, quote="")[1:5, ]

abrupt = read.table('MSIGDB_HUGE_Abruption_Placentae.reformatted.tsv', sep='\t',
                  header=T, quote="")[1:5, ]

all_enrich = rbind(gest, preterm, eclam, abrupt)
all_enrich$phenotype = c(rep('gest', 5), rep('preterm', 5), rep('eclam', 5),
                         rep('abrupt', 5))

all_enrich$Gene_Set_Name = factor(all_enrich$Gene_Set_Name,
  levels=unique(all_enrich$Gene_Set_Name[order(-log10(all_enrich$FDR_q.value))]))

ggplot(all_enrich, aes(x=Gene_Set_Name, y=-log10(FDR_q.value))) + 
  geom_bar(stat='identity', fill='gray', col='black') + coord_flip() + 
  facet_wrap(~phenotype, scales='free', nrow=4) + theme_bw()


draw.quad.venn(50, 339, 449, 517, 15, 20, 17, 66, 65, 244, 7, 6, 15, 45, 5)
dev.off()

draw.quad.venn(100, 100, 100, 100, 26, 18, 14, 32, 33, 73, 12, 9, 14, 28, 9)
dev.off()

# Tottge Fugyre 3


corr_matr = as.matrix(read.table('../ukb_gwas_sumstats/pregnancy_matr_corr.txt', 
                       sep='\t', header=F))

distances = as.dist(1-abs(corr_matr))
hc_res <- hclust(distances)
plot(hc_res)
str(hc_res)

heads = read.table('../ukb_gwas_sumstats/head_line', header=T, sep = ' ')
phen_inds = sapply(strsplit(as.character(colnames(heads[(2:46)*2])), '_'),
                   function(x) x[1])
phen_inds = gsub('\\.', '_', phen_inds)
phen_inds = gsub('X2', '2', phen_inds)
phen_inds = gsub('X4', '4', phen_inds)
phen_inds = gsub('X3', '3', phen_inds)

phenotype_stats = read.table('../ukb_gwas_sumstats/BRB5_GWAS_PREGNANCY/publication/preg_phen_snp_lGC.tsv',
                             sep='\t', header=T)

descrs_ordered = sapply(phen_inds[hc_res$order], 
       function(x) as.character(phenotype_stats
                                [which(phenotype_stats$Phenotype_Code == x), 
                                  'Phenotype_Description'][1]))
phenotype_names = phenotype_stats$Phenotype_Description[]

corrs = melt(corr_matr)
corrs$Var1 = factor(corrs$Var1, levels=unique(corrs$Var1)[hc_res$order])
corrs$Var2 = factor(corrs$Var2, levels=unique(corrs$Var2)[hc_res$order])


# Export as 15 x 11
ggplot(corrs, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_fill_gradientn(colours=jet) + 
  scale_x_discrete(labels=as.character(descrs_ordered)) +
  scale_y_discrete(labels=as.character(descrs_ordered)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


good_traits = c('20002_1221', '4041', 'I9_HYPTENSPREG',
                'PRE_OR_ECLAMPSIA', 'O62', 'O16', 'O14', 'O60', 'O13',
                '20002_1073')

good_traits_indices = which(phenotype_stats$Phenotype_Code %in% good_traits)
good_traits_order = hc_res$order[hc_res$order %in% good_traits_indices]
good_traits_descrs = phenotype_stats$Phenotype_Description[good_traits_order]
good_traits_matrix = corrs[as.character(corrs$Var1) %in% good_traits_indices &
                          gsub('V', '', as.character(corrs$Var2)) %in% good_traits_indices, ]

ggplot(good_traits_matrix, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_fill_gradientn(colours=mypal) + 
  scale_x_discrete(labels=as.character(good_traits_descrs)) +
  scale_y_discrete(labels=as.character(good_traits_descrs)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


#  Plotting lambda vs beritability
phenotype_stats$h2_observed[is.na(phenotype_stats$h2_observed)] = 0.00

pstat = melt(phenotype_stats, id.vars=c('Phenotype_Code'),
             measure.vars=c('lambdaGC', 'h2_observed'))

ggplot(pstat, aes(x=Phenotype_Code, y=value, fill=variable)) + 
  geom_bar(stat='identity', position='dodge', col='black') +
  facet_wrap(~variable, nrow=2) +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=c(mycol1, mycol2))


ggplot(phenotype_stats, aes(x=h2_observed, size=log10_P7, 
                            y=lambdaGC)) + geom_point() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))


# GWAS Catalog histogram

gwcat = read.table('../ukb_gwas_sumstats/BRB5_GWAS_PREGNANCY/publication/gwas_catalog_filtered_pval_05.tsv',
                   sep='\t', header=T)

gwcat$Reported_trait = factor(gwcat$Reported_trait, 
        levels=names(sort(table(as.character(gwcat$Reported_trait)), )))

ggplot(gwcat, aes(x=Reported_trait)) + geom_bar(col='black', fill='gray') +
  coord_flip() + theme_bw()


draw.triple.venn(555, 706, 327, 24, 9, 5, 1)
dev.off()

# For SI tables

regions_data = read.table(
  '../ukb_gwas_sumstats/final_tables/ukb/ukb_genomic_regions_wSNPs.tsv',
                          sep='\t', header=F)

colnames(regions_data) = c('Chromosome', 'Start', 'End', 'Phenotype_Code',
                           'Genes', 'SNPs', 'p_values')

get_descr <- function(phenotypes){
  descrs = c()
  for(phen in phenotypes){
    this_descr = as.character(phenotype_stats[as.character(phenotype_stats$Phenotype_Code) == as.character(phen),
                                              'Phenotype_Description'][1])
    descrs = c(descrs, this_descr)
  }
  return(paste(descrs, sep=',', collapse=""))
}

regions_data$Phenotype_Description = sapply(regions_data$Phenotype_Code,
  function(x) get_descr(strsplit(as.character(x), ',')[[1]]))

write.table(regions_data, 
            file='../ukb_gwas_sumstats/final_tables/SI_Table_4_UKB_loci.tsv', 
            sep='\t', quote=F,col.names=NA)

### Enrichments and curated gene sets graphs

library(dplyr)
library(msigdbr)
library(clusterProfiler)

cp_tab = msigdbr(species='Homo sapiens', category="C2")
cp_tab = cp_tab[cp_tab$gs_subcat != 'CGP', ]
m_t2g <- select(.data = cp_tab, gs_name, gene_symbol)

dg_symbols = read.table('./curated_list/diabetes_gestational.gene_list', heade=F)$V1
dg_enrich = enricher(as.character(dg_symbols), TERM2GENE=m_t2g) %>% as_data_frame()
write.table(dg_enrich, file='CP_MsigDB_curated_Diabetes_Gestational.tsv',
            sep='\t', row.names=F, quote=F)


pe_symbols = read.table('./curated_list/preeclampsia.gene_list', heade=F)$V1
pe_enrich = enricher(as.character(pe_symbols), TERM2GENE=m_t2g) %>% as_data_frame()
write.table(pe_enrich, file='CP_MsigDB_curated_Preeclampsia.tsv',
            sep='\t', row.names=F, quote=F)
            

ptb_symbols = read.table('./curated_list/preterm_birth.gene_list', heade=F)$V1
ptb_enrich = enricher(as.character(ptb_symbols), TERM2GENE=m_t2g) %>% as_data_frame()
write.table(ptb_enrich, file='CP_MsigDB_curated_Preterm_birth.tsv',
            sep='\t', row.names=F, quote=F)
            

pa_symbols = read.table('./curated_list/placental_abruption.gene_list', heade=F)$V1
pa_enrich = enricher(as.character(pa_symbols), TERM2GENE=m_t2g) %>% as_data_frame()
write.table(pa_enrich, file='CP_MsigDB_curated_Placental_Abruption.tsv.tsv',
            sep='\t', row.names=F, quote=F)
            

all_enrich = rbind(dg_enrich[1:5, ], ptb_enrich[1:5, ], 
                   pe_enrich[1:5, ], pa_enrich[1:5, ])
all_enrich$phenotype = c(rep('gest', 5), rep('preterm', 5), rep('eclam', 5),
                         rep('abrupt', 5))

all_enrich$ID = factor(all_enrich$ID, 
                       levels=unique(all_enrich$ID[order(-log10(all_enrich$p.adjust))]))

ggplot(all_enrich, aes(x=ID, y=-log10(p.adjust))) + 
  geom_bar(stat='identity', fill='gray', col='black') + coord_flip() + 
  facet_wrap(~phenotype, scales='free', nrow=4) + theme_bw()

dev.off()

draw.quad.venn(42, 170, 284, 233, 9, 14, 12, 41, 34, 103, 6, 4, 9, 24, 4)
dev.off()

draw.quad.venn(81, 181, 292, 336, 40, 53, 52, 92, 93, 214, 33, 32, 45, 75, 30)
dev.off()

# Uniques

all_all_enrich = rbind(dg_enrich, pa_enrich, pe_enrich, ptb_enrich)

unique_genes = read.table('unique_gene_list.txt', header=F)

unique_pathways = list()
for (i in unique(unique(as.character(unique_genes$V2)))) {
  print(i)
  this_enrich = enricher(as.character(unique_genes[unique_genes$V2 == i, 'V1']), 
                         TERM2GENE=m_t2g) %>% as_data_frame()
  write.table(this_enrich, file=paste0('CP_unique_genes_', i, '.tsv'),
              sep='\t', quote=F, row.names=F)
  uniq_path = this_enrich$ID[!(this_enrich$ID %in% all_all_enrich$ID)]
  if (length(uniq_path) == 0) {
    uniq_path = c('NULL')
  }
  unique_pathways[i] = uniq_path
}
