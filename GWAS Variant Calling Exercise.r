The task is to retrieve GWAS summary statistics from following datasets : 
		-Obesity (Bradfield JP Nature Genetics 2012)
		-Insulin measurement (Strawbridge RJ Diabetes 2011)
		-Fasting blood glucose measurement (Dupuis J Nature Genetics 2010)
For each summary stat take only the variants with:
		-P-value lower than 0.01
		-Allele Frequency greater than 1%
For variants from each summary dataset, retrieve following additional information:
		-Chromosome
		-Position on GRCh37
		-Non Effective allele
		-Effective allele
		-ALT allele Frequency
		-Beta
		-Beta stderr
		-P-value
The second required output is a CSV file with a list of variants that are present in all three summary stats.



#Retrieving the summary statisitcs for Obesity, Insulin measurement, Fasting blood glucose measurement
gwahits_Obs = data.frame(read.table('EGG_Obesity_Meta_Analysis_1.txt',header=TRUE,sep=' ',as.is=TRUE))

gwahits_FastIn = data.frame(read.table('MAGIC_ln_FastingInsulin.txt',header=TRUE,sep='\t',as.is=TRUE))
gwahits_proIn = data.frame(read.table('MAGIC_proinsulin_for_release_HMrel27.txt',header=TRUE,sep='\t',as.is=TRUE))

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

#gwahits_FastInMafPval<-gwahits_FastIn[gwahits_FastIn$maf>0.01]


#Selecting data with selection critirea maf>0.01 and pvalue<0.01
#The Obesity dataset didnt contain maf value
library(dplyr)
gwahits_FastInMafPval<-filter(gwahits_FastIn, maf > 0.01 & pvalue<0.01)
gwahits_ProInMafPval<-filter(gwahits_proIn, maf > 0.01 & pvalue<0.01)





if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantTools")
BiocManager::install("biomaRt")
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
#######################
#Retrieving from biomaRt
library(biomaRt)
# biomaRt_2.30.0, R version 3.3.2 (2016-10-31)
#getBM function has a timeout so queries are broken down to sets of 1/10
sn<-gwahits_ProInMafPval[24000:27000,]
sn$snp
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
snp.db <- useMart("ENSEMBL_MART_SNP",host="grch37.ensembl.org", dataset="hsapiens_snp")
listAttributes(snp.db)
nt.biomart_Pro9 <- getBM(c("refsnp_id","allele","chr_name","chrom_start",                   
                      "chrom_strand",
                      "allele_1","minor_allele","minor_allele_freq"),
                    filters="snp_filter",
                    values=sn$snp,
                    mart=snp.db)


###############Merging the fetched results




nt.biomart_Pro_Final <- rbind(nt.biomart_Pro1,nt.biomart_Pro2,nt.biomart_Pro3,nt.biomart_Pro4,nt.biomart_Pro5,nt.biomart_Pro6,
             nt.biomart_Pro7,nt.biomart_Pro8,nt.biomart_Pro9)

nt.biomart_Fast_Final <- rbind(nt.biomart_Fast1,nt.biomart_Fast2,nt.biomart_Fast3,nt.biomart_Fast4,nt.biomart_Fast5,nt.biomart_Fast6,
                              nt.biomart_Fast7,nt.biomart_Fast8,nt.biomart_Fast9,nt.biomart_Fast10)

write.csv(nt.biomart_Pro_Final,'nt.biomart_Pro_Final.csv')
write.csv(nt.biomart_Fast_Final,'nt.biomart_Fast_Final.csv')

Pro_Combined<-merge(x = gwahits_ProInMafPval, y = nt.biomart_Pro_Final,by.x = "snp", by.y = "refsnp_id", all = TRUE)
Fast_Combined<-merge(x = gwahits_FastInMafPval, y = nt.biomart_Fast_Final,by.x = "snp", by.y = "refsnp_id", all= TRUE)



#Saving the dataset
write.csv(Pro_Combined,'Insulin_Measurement_Pro_Results.csv')
write.csv(Fast_Combined,'Fast_blood_measurement_results.csv')

library(dplyr)
Pro_Fast_Combined<-merge(x = Fast_Combined, y = Pro_Combined,by.x= "snp",by.y="snp", all = TRUE)


#Inner join the snp data from tables from both papers
Pro_Fast_Combined1<-inner_join(Fast_Combined, Pro_Combined, by = "snp")


write.csv(Pro_Fast_Combined1,'Pro_Fast_Combined1.csv')


#########################Preprocessing on Obesity dataset to add maf ###################################################





gwahits_ObsMafPval<-filter(gwahits_Obs, PVAL<0.01)
sn<-gwahits_ObsMafPval[24000:28579,]#Make subset of the gwahits_ObsMafPval 
sn$SNP

snp.db <- useMart("ENSEMBL_MART_SNP",host="grch37.ensembl.org", dataset="hsapiens_snp")
nt.biomart_Obs9 <- getBM(c("refsnp_id","allele","chr_name","chrom_start",                   
                           "chrom_strand",
                           "allele_1","minor_allele","minor_allele_freq","p_value"),
                         filters="snp_filter",
                         values=sn$SNP,
                         mart=snp.db)

nt.biomart_Obs_Final <- rbind(nt.biomart_Obs1,nt.biomart_Obs2,nt.biomart_Obs3,nt.biomart_Obs4,nt.biomart_Obs5,nt.biomart_Obs6,
                              nt.biomart_Obs7,nt.biomart_Obs8,nt.biomart_Obs9)




nt.biomart_Obs_Final

#Checking Maf criterian after fetching MAF data
nt.biomart_Obs_Final_with_Maf<-filter(nt.biomart_Obs_Final, minor_allele_freq>0.01)
nt.biomart_Obs_Final_with_Maf
nt.biomart_Obs_Final
Obs_Combined1<-merge(x = gwahits_ObsMafPval, y = nt.biomart_Obs_Final_with_Maf,by.x = "SNP", by.y = "refsnp_id", all.x= TRUE)

write.csv(Obs_Combined1,'Obesity_Results.csv')

################################################Combining for Common SNPs####################



Obs_Combined2<-Obs_Combined1
#changing column RsID name to match tables
names(Obs_Combined2)[1] <- "snp"
#Combining with Obesity, ProInsulin and Fating dataset
All_combined<-inner_join(Obs_Combined2, Pro_Fast_Combined1, by = "snp")
All_combined_results<-All_combined[,1:13]
All_combined_results<-na.omit(All_combined_results)
write.csv(All_combined,'All_common_snps.csv')

library(tidyverse)
All_combined_results_table<-distinct(All_combined_results, .keep_all = TRUE)
write.csv(All_combined_results_table,'All_Stats_common_snp_results_table.csv')









###########################END#######################################
#####################                      ###############################
#####################                      ###########################
