BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
library(org.Hs.eg.db)

# convert ENSEMBL to SYMBOL using AnnotationDbi package
table_hs = AnnotationDbi::select(
  org.Hs.eg.db, 
  keytype="ENSEMBL",
  keys=eng_vec, 
  columns=c("SYMBOL") 
)

# trying to fill NA values with values from ensg_convert.csv which is a list of genes in ENG and their 
# symbol. i got this list from https://www.biotools.fr/mouse/ensembl_symbol_converter .
# all_genes is the final couplind table.

not_NA<-table_hs[!is.na(table_hs$SYMBOL),]
not_NA<-not_NA[!duplicated(not_NA$SYMBOL),]

NA_genes<-table_hs[is.na(table_hs$SYMBOL),]
ensg_convert_csv <- read.csv(file = 'ensg_convert.csv',check.names=FALSE)
NA_genes<-merge(NA_genes,ensg_convert_csv, by='ENSEMBL',all=TRUE)
NA_genes<-NA_genes[ , c("ENSEMBL","symbol")]
colnames(NA_genes)[2] <- "SYMBOL"

blank<-NA_genes[(NA_genes$SYMBOL==""),]
NA_genes<-NA_genes[!(NA_genes$SYMBOL==""),]
NA_genes<-NA_genes[!duplicated(NA_genes$SYMBOL),]
all_genes<-rbind(not_NA, NA_genes)
all_genes<-all_genes[!duplicated(all_genes$SYMBOL),]


# The conversion 
colnames(all_genes)[1] <- "ENSEMBL"
A_ENG_data<-merge(all_genes,ENG_data, by='ENSEMBL',all=TRUE)
A_ENG_data$ENSEMBL<- NULL

rm(table_hs,NA_genes,blank,ensg_convert_csv, not_NA)

