# Какие-то полезные ссылки и слова,ну по которым сходу можно добраться до осложнений беременности... 
# Maternal Hypothyroidism
# гестационный диабет, эклампсия
# preterm labor
# placental abruption
# 
# <https://icd.codes/icd10cm/chapter15>
# Code Range 	Section Description 
# O00-O08 	Pregnancy with abortive outcome
# O09 	Supervision of high risk pregnancy
# O10-O16 	Edema, proteinuria and hypertensive disorders in pregnancy, childbirth and the puerperium
# O20-O29 	Other maternal disorders predominantly related to pregnancy
# O30-O48 	Maternal care related to the fetus and amniotic cavity and possible delivery problems
# O60-O77 	Complications of labor and delivery
# O80-O82 	Encounter for delivery
# O85-O92 	Complications predominantly related to the puerperium
# O-O9A 	Other obstetric conditions, not elsewhere classified
# 
# 
# 
# 20002_1559 - self-reported miscarriage
# 3839 - number of miscarriages
# 2774 -Ever had stillbirth, spontaneous miscarriage or termination
# 3829 - Number of stillbirths
# Они все довольно нулевые в плане наследуемости по оценкам на ukbb
# Но добавить в датасет их я думаю, можно
# Я это тебе говорю потому, что как раз глотовцев в первую очередь интересуют как раз miscarriage и termination в этих терминах
# То есть поскольку они хотят писать этот полуобзор-полуанализ открытых данных

# rename "s/bgz/gz/" *.bgz # one-liner to make .bgz files readable for data.table::fread

setwd("~/Documents/Bioinf/BRB5/")
library(dplyr)
library(readxl)
library(stringr)
library(tidyr) 
library(data.table)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(qqman)
library(GenABEL) 
library(CMplot)
library(viridis)
library(pbapply)  
library(phenoscanner)
library(eulerr)

# I had a problem with an installation of "GenABel" package (R version 3.6.1 (2019-07-05))
# So this is the way to fix it
# install.packages("devtools")
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz")
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz")


uk_biobank <- read_excel("UKBB GWAS Imputed v3 - File Manifest Release 20180731.xlsx",
                         sheet = 2,na = "N/A",col_names = T)
# выкидываю все что без imputed.v3
uk_biobank_v3 <- subset(uk_biobank,grepl(uk_biobank$File,pattern="imputed.v3") & uk_biobank$Sex == "female") 

# покопался в данных uk biobank нашел какой-то файл phenotypes.female.tsv
# там вроде есть коды фенотипов по которым я пробовал сабсетить по ключевым словам 
phenotype_fem <- read.delim("phenotypes.female.tsv")

# # что то связанное с осложнениями при беременности, МОЖНО ДОБАВИТЬ БОЛЬШЕ ДАННЫХ!
# pregn_prob <- subset(phenotype_fem,grepl(description,pattern = "pregn") | grepl(description,pattern = "labor") |
#                          grepl(description,pattern = "ypothyroidism") | grepl(description,pattern = "clamps") | 
#                          grepl(description,pattern = "estational"))
# 
# после добавлял в ручную некоторые отдельные коды и в итоге wget'ом выкачал все данные -> my_data


preg_phen <- read.table("data/my_data.txt",header = F)
preg_phen$`Phenotype Code` <- gsub(preg_phen$V9,pattern = "\\..*",replacement="")
colnames(preg_phen)[9] <- "file"
preg_phen <- plyr::join(preg_phen[,-c(1:8)], uk_biobank_v3[,c(1,2,6)],type="inner")


rm(uk_biobank,uk_biobank_v3,phenotype_fem)

"For each dataset it calculates:
lambdaGC, number of all snp, high confident snp and snp with pval < 1e-06,1e-07,1e-08
but it put them to one column, so later we'll separate it to multiple columns"

preg_phen$result <- pbsapply(preg_phen$file, function(x){
  
  data <- fread(paste("data/",x, sep=""))
 
  data %>% drop_na() %>% count() -> raw_snp1
  
  data[low_confidence_variant==FALSE] %>% drop_na() -> data
  data %>% count() -> highconf_snp
  
  data %>%  
    dplyr::select(pval) %>%
    summarise(lambdaGC = estlambda(data = pval,plot=F,method = "median")[1]) -> lambdaGC
    as.vector(unlist(lambdaGC)) %>% round(.,digits = 3) -> lambdaGC

  data[pval < 1e-08] %>% count() -> neglog10_P8
  data[pval < 1e-07] %>% count() -> neglog10_P7
  data[pval < 1e-06] %>% count() -> neglog10_P6
  
  total <- paste(raw_snp1,highconf_snp,lambdaGC,neglog10_P6,neglog10_P7,neglog10_P8, sep = ";")
  return(total)
})

# It creates several columns to put there the corresponding information from 'result' column
# and deletes result column in the end
preg_phen$raw_snp <- NA
preg_phen$highconf_snp <- NA
preg_phen$lambdaGC <- NA
preg_phen$log10_P8 <- NA
preg_phen$log10_P7 <- NA
preg_phen$log10_P6 <- NA

preg_phen[,c(6:ncol(preg_phen))] <- data.frame(str_split_fixed(preg_phen$result,pattern = ';',n=Inf))
preg_phen %>% dplyr::select(-result) -> preg_phen
colnames(preg_phen) <- gsub(colnames(preg_phen),pattern = " ",replacement = "_")
# this dataframe contains the data (phenotypes)
# I used for counting SNP with different pval thresholds (-log(p) = [6,7,8]), 
# lambdaGC, raw and high confident SNP numbers
fwrite(preg_phen,"preg_phen_snp_lGC.tsv",sep = "\t") 

#
preg_phen <- read.table("preg_phen_snp_lGC.tsv",header = T,sep = "\t")

# Visualization of preg_phen (pval < 1e-07)
plot1 <- ggplot(preg_phen, aes(Phenotype_Code,lambdaGC))+
  geom_point(aes(size=as.numeric(preg_phen$log10_P7),description = Phenotype_Description, snp = log10_P7),col=brewer.pal(7,"Dark2")[1],alpha=0.9)+
  ylab("Genomic inflation factor")+
  xlab("Phenotype code")+
  labs(size="Number of significant \nSNP (1-e07)")+
  geom_hline(yintercept = 1)+ylim(limits = c(0.9,1.2))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45,size = 10,face = "plain"))+
  scale_size_continuous(range = c(1.5,8))+ylim(NA,1.05)
plot1
# ggplotly for more convinient analysis (to see Phenotype_Code, Phenotype_Description and number of snp with the correspoding threshold)
ggplotly(plot1,tooltip = c("Phenotype_Code","description","snp"))

# Visualization of preg_phen (pval < 1e-08)
plot2 <- ggplot(preg_phen, aes(Phenotype_Code,lambdaGC))+
  geom_point(aes(size=as.numeric(preg_phen$log10_P8),description = Phenotype_Description, snp = log10_P8),col=brewer.pal(7,"Dark2")[3],alpha=0.9)+
  ylab("Genomic inflation factor")+
  xlab("Phenotype code")+
  labs(size="Number of significant \nSNP (1-e08)")+
  geom_hline(yintercept = 1)+ylim(limits = c(0.9,1.2))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45,size = 10,face = "plain"))+
  scale_size_continuous(range = c(1.5,8))+ylim(NA,1.05)
plot2
ggplotly(plot2,tooltip = c("Phenotype Code","description","snp"))

# we see, that only 4 datasets have a signal (number of SNP)
# now we're going to analyse them more properly (QQ and Manhattan plots)

# Function for edition datasets for Manhattan and QQ plots
gwas_edit <- function(data){
  data[low_confidence_variant==F][,c(1,12)] %>% drop_na()-> data
  colnames(data)[1] <- "SNP"
  data$Chromosome <- NA
  data$Position <- NA
  data[,c(3,4)] <- as.data.frame(str_split_fixed(data$SNP,":",Inf)[,c(1,2)])
  data <- data[,c(1,3,4,2)]
  data$Position <- sapply(data$Position, as.numeric)
  return(data)
}
# Function for manhattan plot 
# WARNING! It creates Manhattan plot from pval > 0.01 to save time
manhattan_plot <- function(data,pval_threshold,outfile){
  CMplot((data %>% filter(pval < pval_threshold)),
         plot.type = "m",band = 0,LOG10 = TRUE,
         threshold = c(1e-6,1e-7,1e-8),threshold.col = c("blue","green","red"),
         amplify = F,ylim = c(2,9),cex = 0.8,
         file="jpg",memo="",dpi=300,file.output=outfile,verbose=TRUE,width=7,height=6)
}
# Function for QQ plot
qq_plot <- function(data,outfile){
  CMplot(data,plot.type="q",conf.int=TRUE,box=FALSE,file="jpg",memo="",dpi=300,
         ,file.output=outfile,verbose=TRUE,width=5,height=5)
}



# Processing O46 dataset (creating all plots and subsetting for pval < 1e-05)
O46 <- fread("data/O46.gwas.imputed_v3.female.tsv.gz")
O46_p <- gwas_edit(O46) # dataframe for plotting (Manhattan or QQ)
subset(O46_p, O46_p$pval <= 1e-5) -> snp_O46
manhattan_plot(O46_p,0.01,T)
colnames(O46_p)[4] <- "O46 phenotype"
qq_plot(O46_p,outfile = T)
rm(O46,O46_p)

# Processing O26 dataset (creating all plots and subsetting for pval < 1e-05)
O26 <- fread("data/O26.gwas.imputed_v3.female.tsv.gz")
O26_p <- gwas_edit(O26) # dataframe for plotting (Manhattan or QQ)
subset(O26_p, O26_p$pval <= 1e-5) -> snp_O26
manhattan_plot(O26_p,0.01,T)
colnames(O26_p)[4] <- "O26 phenotype"
qq_plot(O26_p,outfile = T)
rm(O26,O26_p)

# Processing O69 dataset (creating all plots and subsetting for pval < 1e-05)
O69 <- fread("data/O69.gwas.imputed_v3.female.tsv.gz")
O69_p <- gwas_edit(O69) 
subset(O69_p, O69_p$pval <= 1e-5) -> snp_O69
manhattan_plot(O69_p,0.01,T)
colnames(O69_p)[4] <- "O69 phenotype"
qq_plot(O69_p,outfile = T)
rm(O69,O69_p)

# Processing I9_HYPTENSPREG dataset (creating all plots and subsetting for pval < 1e-05)
I9_HYPTENSPREG <- fread("data/I9_HYPTENSPREG.gwas.imputed_v3.female.tsv.gz")
I9_HYPTENSPREG_p <- gwas_edit(I9_HYPTENSPREG) 
subset(I9_HYPTENSPREG_p, I9_HYPTENSPREG_p$pval <= 1e-5) -> snp_I9_HYPTENSPREG
manhattan_plot(I9_HYPTENSPREG_p,0.01,T)
colnames(I9_HYPTENSPREG_p)[4] <- "I9_HYPTENSPREG phenotype"
qq_plot(I9_HYPTENSPREG_p,outfile = T)
rm(I9_HYPTENSPREG,I9_HYPTENSPREG_p)

# creating summary table for these 4 datasets
snp_I9_HYPTENSPREG$Dataset <- "I9_HYPTENSPREG" 
snp_O26$Dataset <- "O26" 
snp_O46$Dataset <- "O46" 
snp_O69$Dataset <- "O69" 

ukbiobank_summ <- rbind(snp_I9_HYPTENSPREG,snp_O26,snp_O46,snp_O69)
write.csv(ukbiobank_summ[,c(1,4)],"snp_ukbiobank_summ.csv")





# Processing O26 dataset (creating all plots and subsetting for pval < 1e-05)
O46 <- fread("data/O46.gwas.imputed_v3.female.tsv.gz") 
O46 %>% filter(low_confidence_variant==F & O46$pval <= 1e-05) -> snp_O46
rm(O46)
# Processing O46 dataset (creating all plots and subsetting for pval < 1e-05)
O26 <- fread("data/O26.gwas.imputed_v3.female.tsv.gz") 
O26 %>% filter(low_confidence_variant==F & O26$pval <= 1e-05) -> snp_O26
rm(O26)
# Processing I9_HYPTENSPREG dataset (creating all plots and subsetting for pval < 1e-05)
I9_HYPTENSPREG <- fread("data/I9_HYPTENSPREG.gwas.imputed_v3.female.tsv.gz") 
I9_HYPTENSPREG %>% filter(low_confidence_variant==F & I9_HYPTENSPREG$pval <= 1e-05) -> snp_I9_HYPTENSPREG
rm(I9_HYPTENSPREG)
# Processing O69 dataset (creating all plots and subsetting for pval < 1e-05)
O69 <- fread("data/O69.gwas.imputed_v3.female.tsv.gz")
O69 %>% filter(low_confidence_variant==F & O69$pval <= 1e-05) -> snp_O69
rm(O69)
# creating a column with dataset information for each variant
snp_I9_HYPTENSPREG$Dataset <- "I9_HYPTENSPREG" 
snp_O26$Dataset <- "O26" 
snp_O46$Dataset <- "O46" 
snp_O69$Dataset <- "O69" 

total <- rbind(snp_O46,snp_O26,snp_I9_HYPTENSPREG,snp_O69)
write.csv(total,"total_ukbiobank_pregn.csv")

#
total <- read.csv("total_ukbiobank_pregn.csv")

variants <- fread("variants.tsv.gz",nrows = 50) # dataset from UKB with metadata for each SNP
plyr::join(total,variants[,c(1:9)],type="inner") -> total
rm(variants)

# saving data for LSEA analysis for all 4 phenotypes from UKB
colnames(total)[c(14,15,18,16,17,12)] <- c("CHR","COORDINATE","RSID","REF","ALT","PVAL")
total[,c("CHR","COORDINATE","RSID","REF","ALT","PVAL","Dataset")] -> ukb_total_lsea

# filtering snp wich don't have rd-id
ukb_total_lsea %>% filter(grepl(.$RSID,pattern='rs')) -> ukb_total_lsea_filt

# exporting tables for 4 UKB phenotypes separately for LSEA
sapply(unique(ukb_total_lsea_filt$Dataset),function(x){
  data <- subset(ukb_total_lsea_filt,ukb_total_lsea_filt$Dataset==x)
  full_name = paste0(noquote(x),".csv",collapse = "")
  write.csv(data[,c("CHR","COORDINATE","REF","ALT","PVAL")],full_name,row.names = F)
})

# exporting tables for all phenotypes separately for LSEA
ukb_total_lsea_filt %>% dplyr::select(-Dataset) %>% write.csv(.,'UKB_total.csv',row.names = F)

###################################################################################################################################
###################################################################################################################################
#######       ###    ###     #####    ##      ######         #########       #####      ###           #############################
#######   ########   ##   #   ###   ###   ##   #####   ###############   ########   ##   ######   #################################
#######   ##   ####      ###   #   ###   ####   ####         #########   #######   ####   #####   #################################
#######   ###  #####    #####     ###            #########   #########   ######            ####   #################################
#######        ######  #######   ###    #####    ###         #########       #   ########   ###   #################################
###################################################################################################################################
###################################################################################################################################

# Now we're working with GWAS Catalog 
# gwas_catalog table was manually obtained from GWAS Catalog website
# There we just took all snps which somehow corresponded to pregnancy problems
gwas_catalog <- read_excel("gwas_catalog_rawtable.xlsx",col_names = T)
colnames(gwas_catalog) <- gsub(colnames(gwas_catalog),pattern=" ",replacement = "_")
# Calculation unique and total snp in the whole table
gwas_catalog %>% summarise(total_snp=nrow(.),unique_snp=length(unique(gwas_catalog$Variant_and_risk_allele)))

# this is a file with selected phenotypes from the GWAS Catalog more precisely, removing excess phenotypes
gwas_traits <- read.csv("filt.csv")
gwas_traits$V3 <- paste(gwas_traits$V1,gwas_traits$V2,sep = "_")
gwas_catalog <- subset(gwas_catalog,gwas_catalog$Reported_trait %in% gwas_traits$V1 & gwas_catalog$`Trait(s)` %in% gwas_traits$V2)
gwas_catalog <- subset(gwas_catalog, gwas_catalog$Reported_trait %in% gwas_traits$V1 & gwas_catalog$`Trait(s)` %in% gwas_traits$V2)
gwas_traits <- as.data.frame(table(gwas_catalog$Reported_trait))

# editing of gwas_catalog dataframe 
gwas_catalog$RSID <- as.data.frame(str_split_fixed(gwas_catalog$Variant_and_risk_allele, pattern = "-",n=Inf))[,1]
colnames(gwas_catalog)[2] <- "PVAL"
gwas_catalog$REF <- as.data.frame(str_split_fixed(as.data.frame(str_split_fixed(gwas_catalog$snp_info,pattern=";",n=Inf))[,2],pattern="\\/",Inf))[,1]

# calculation of phenotypes from the GWAS Catalog and the number of SNPs which they account for
ggplot(gwas_traits,aes(reorder(Var1,Freq),Freq))+
  geom_bar(stat = 'identity',fill=brewer.pal(3,'Set1')[2])+
  coord_flip()+
  theme_bw()+ggtitle("GWAS Catalog traits\nrelated to pregnancy problems")+
  ylab("Counts")+xlab("GWAS Catalog phenotypes")+
  theme(axis.text.y = element_text(size=14,colour="black"),
        title = element_text(size=20),axis.text.x=element_text(size=14,colour='black'))


# OUR COLABORATORS ASKED FOR THIS TABLE TO MAKE SOME COMBINATIONS BETWEEN THESE PHENOTYPES
write.csv(gwas_traits,"GWAS_for_collab.csv")
getwd()


# correction PVAL column (making it numeric type)
d <- as.data.frame(str_split_fixed(gwas_catalog$PVAL,pattern=" x ",Inf))
d$V1 <- as.numeric(as.character(d$V1))
d <- cbind(d[,1],as.data.frame(str_split_fixed(d$V2,pattern="-",Inf))) 
d[,c(2,3)] <- apply(d[,c(2,3)],c(1,2),as.character) 
d[,c(2,3)] <- apply(d[,c(2,3)],c(1,2),as.numeric) 
d$V1 <- 10**(-d$V2)
gwas_catalog$PVAL <- d$V1
rm(d)

# PHENOSCANNER PACKAGE R SNP annotation 
gwas_catalog$NUC_EXCHANGE <- NA
gwas_catalog$REGION <- NA

gwas_catalog %>%
  dplyr::mutate(NUC_EXCHANGE = pbsapply(.$RSID, function(x){
    res <- phenoscanner(snpquery=x,catalogue = "GWAS")
    if (!(grepl(x,pattern = 'chr')) | all(dim(res$snps)!=0)){
      return(paste0(res$snps[,c(8,9)],collapse="/"))
    } else {
      return("no_data")
    }
  }),
  REGION = pbsapply(.$RSID, function(x){
    res <- phenoscanner(snpquery=x,catalogue = "GWAS")
    if (all(dim(res$snps)!=0)){
      return(paste0(res$snps$consequence,collapse=";"))
    } else {
      return("no_data")
    }
  })
  ) -> gwas_catalog


# counting SNPs with manual annotation manual annotation (25 in total)
gwas_catalog %>% 
  filter(grepl(.$NUC_EXCHANGE,pattern='-')) %>% count()

# counting failed SNP id which contain anything but rsID  (14 in total)
gwas_catalog %>% 
  filter(grepl(.$Variant_and_risk_allele,pattern='rs')==F) %>% count()

# exporting gwas_catalog to be safe
write.csv(gwas_catalog,"gwas_catalog_annotated.csv")
# then we manually check the SNPs that have crap in REF/ALT columns (chr etc)
# adding REF and ALT from GWAS Catalog

# Now preprocessing of variants from gwas catalog for LSEA
gwas_catalog <- read_excel("gwas_catalog_annotated.xlsx")
gwas_catalog$Location <- gsub(gwas_catalog$Location,pattern="X",replacement = "23") # заменяю Х хромосому на 23
gwas_catalog$REF <- gwas_catalog$REF_man_corr
gwas_catalog$REF <- ifelse(is.na(gwas_catalog$REF),
                                    as.character(as.data.frame(str_split_fixed(gwas_catalog$NUC_EXCHANGE,pattern = "/",Inf))[,1]),
                                    gwas_catalog$REF)
gwas_catalog$ALT <- gwas_catalog$ALT_man_corr
gwas_catalog$ALT <- ifelse(is.na(gwas_catalog$ALT) & gwas_catalog$NUC_EXCHANGE!='no_data' ,
                                    as.character(as.data.frame(str_split_fixed(gwas_catalog$NUC_EXCHANGE,pattern = "/",Inf))[,2]),
                                    gwas_catalog$ALT)
gwas_catalog$ALT <- ifelse(grepl(gwas_catalog$ALT,pattern='/'),
                                    as.character(as.data.frame(str_split_fixed(gwas_catalog$ALT,pattern = "/",Inf))[,2]),
                                    gwas_catalog$ALT)
gwas_catalog %>% 
  summarise(non_valid_snp = subset(.,grepl(.$REF,pattern='-') | grepl(.$REF,pattern='no_data') | is.na(.$REF) | 
                                 grepl(.$ALT,pattern='-') | grepl(.$ALT,pattern='no_data') | is.na(.$ALT) |
                                   is.na(.$COORDINATE)) %>% nrow(.), 
            valid_snp = nrow(.)-non_valid_snp,
            total_snp = nrow(.))

gwas_catalog %>%
  filter(.,grepl(.$REF,pattern='-')==F & grepl(.$REF,pattern='no_data')==F & is.na(.$REF)==F &
           grepl(.$ALT,pattern='-')==F & grepl(.$ALT,pattern='no_data')==F & is.na(.$ALT)==F & 
           is.na(.$COORDINATE)==F) -> gwas_catalog_filtered

gwas_catalog_filtered$CHR <- as.data.frame(str_split_fixed(gwas_catalog_filtered$Location,pattern = ":",Inf))[,1]

# Export gwas catalog snp for LSEA (total and separated by 25 selected phenotypes)
gwas_catalog_filtered %>% dplyr::select(c(CHR,COORDINATE,RSID,REF,ALT,PVAL)) %>% unique(.) -> gwas_catalog_filtered_lsea

# Total snp from GWAS Catalog for LSEA
write.csv(gwas_catalog_filtered_lsea,'~/Documents/Bioinf/BRB5/RESULTS/gwas_catalog_filtered_lsea.csv',row.names = F)

setwd("~/Documents/Bioinf/BRB5/RESULTS/tables_gwascat/")
# separated phenotypes for LSEA
sapply(gwas_traits$Var1,function(x){
  data <- subset(gwas_catalog_filtered,gwas_catalog_filtered$Reported_trait==x)
  # removing 'bad' symbols for file names
  left_name <- gsub(x,pattern='\\(',replacement = "")
  left_name <- gsub(left_name,pattern='\\)',replacement = "")
  left_name <- gsub(left_name,pattern='-',replacement = "_")
  left_name <- gsub(left_name,pattern=' ',replacement = "_")
  left_name <- gsub(left_name,pattern='\\/',replacement = "_")
  full_name = paste0(noquote(left_name),".csv",collapse = "")
  write.csv(data[,c("CHR","COORDINATE","RSID","REF","ALT","PVAL")],full_name,row.names = F)
})

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# UK Biobank and GWAS comparison
plot(venn(list("UK Biobank" = unique(ukb_total_lsea_filt$RSID),
                       "GWAS Catalog" = unique(gwas_catalog_filtered$RSID))),
             main = "UK Biobank            GWAS Catalog",
             fill = brewer.pal(3,"Set1"), alpha = 0.5, 
             edges = T, 
             lty = 1,
             lwd = 1.5, 
             labels = F, 
             legend = list(cex = 1.5), 
             quantities = list(cex = 2.5,font=1)) 