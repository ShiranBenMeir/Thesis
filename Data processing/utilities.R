
# This function merges the raw FPKM data with the clinical data and metadata (more id's) of the patients and ruturns the merged data.
MERGED_DATA_generation<- function(data){
  gene_data<- data
  rm(data)
  raw_metadata <- read.csv(file = 'metadata.csv')
  raw_clinical <- read.csv(file = 'clinical_data.csv')
  raw_clinical$X<-NULL
  
  #remove "_fpkm from submitter_id column"
  submitter_id=7 
  for (j in 1:nrow(raw_metadata)){
    splitedVar <- strsplit(raw_metadata[j,submitter_id], "_fpkm")[[1]]
    raw_metadata[j,submitter_id]<-splitedVar
  }
  
  
  raw_clinical<-raw_clinical[!duplicated(raw_clinical[,1]),]
  raw_metadata<-raw_metadata[!duplicated(raw_metadata[,8]),]
  
  MERGED_DATA<-merge(raw_metadata, raw_clinical,by='case_id', all.dataframe_1=TRUE)
  MERGED_DATA<- na.omit(MERGED_DATA)   #remove NA
  MERGED_DATA <- MERGED_DATA[order(rownames(MERGED_DATA)),]
  
  return(MERGED_DATA)
}

# This function chooses specific wanted columns- submitter_id, project_id,days_to_death, vital_status, year_of_diagnosis.
# The function also orgenaizes the data and return it as- TCGA_DATA.
# The clinical data is taken from TCGA
generate_TCGA_DATA_from_TCGA_clinical_data<- function(MERGED_DATA,data){
  selected_columns<- dplyr::select(MERGED_DATA, submitter_id, project_id,days_to_death, vital_status, year_of_diagnosis)
  selected_columns <- data.frame(selected_columns[,-1], row.names = selected_columns[,1],check.names=FALSE)
  TCGA_DATA<-merge(selected_columns,data, by='row.names',all=TRUE)
  TCGA_DATA <- data.frame(TCGA_DATA[,-1], row.names = TCGA_DATA[,1],check.names=FALSE)
  TCGA_DATA<- na.omit(TCGA_DATA)
  return(TCGA_DATA)
}


# This function creates TCGA_DATA using Ceccarelli clinical data. the funtion returns TCGA_DATA.
generate_TCGA_DATA_from_Ceccarelli <- function(MERGED_DATA,data){
  IDH_data <- read.csv(file = 'NIHMS746836-supplement-8-first chart.csv',   header= TRUE)
  colnames(IDH_data)[1] <- "case_submitter_id"
  MERGED_IDH_clinical<-merge(MERGED_DATA, IDH_data,by='case_submitter_id', all.dataframe_1=TRUE)
  selected_columns<- dplyr::select(MERGED_IDH_clinical, submitter_id,
                                   "Survival..months.",
                                   "Transcriptome.Subtype",
                                   "IDH.codel.subtype",
                                   "IDH.status",
                                   "ESTIMATE.immune.score",
                                   "gender",
                                   Grade,
                                   project_id, 
                                   "Vital.status..1.dead.",
                                   "Histology",
                                   "ESTIMATE.stromal.score")
  
  selected_columns <- data.frame(selected_columns[ ,-1], row.names = selected_columns[,1],check.names=FALSE)
  TCGA_DATA<-merge(selected_columns,data, by='row.names',all=TRUE)
  TCGA_DATA <- data.frame(TCGA_DATA[,-1], row.names = TCGA_DATA[,1],check.names=FALSE)
  TCGA_DATA<- na.omit(TCGA_DATA)  
}

# This function uses RNAseq results from chaya and extract the significant genes that are A1/A2 associated.
# The function returns a list of: (1) vector of significant A1/A2 genes (2) a df of all significant genes with
# their A1/A2 label.
RNA_seq_data_processing<- function(){
  library(data.table)
  
  astrocytic_df2 <- read.csv(file = 'Human_astrocytic_genes.csv',check.names=FALSE)
  
  astrocytic_df<-astrocytic_df2
  astrocytic_df<-astrocytic_df[(astrocytic_df$PValue<0.05),]
  astrocytic_df<-astrocytic_df[(abs(as.numeric(astrocytic_df$log2FC))>1),]
  
  astrocytic_df <- data.frame(astrocytic_df[,-1], row.names = astrocytic_df[,1],check.names=FALSE) 
  
  #change names in astrocytic_data (not in TCGA)
  rownames(astrocytic_df)[rownames(astrocytic_df) == "MRVI1"] <- "IRAG1"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "H2AFY2"] <- "MACROH2A2"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "LYPLAL1-AS1"] <- "AL513283.1"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "C8ORF4"] <- "TCIM"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "KIAA1644"] <- "SHISAL1"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "MARCH3"] <- "MARCHF3"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "LINC00982"] <- "PRDM16-DT"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "LINC02170"] <- "ARLNC1"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "C11orf44"] <- "LINC02873"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "LINC01658"] <- "FAM230I"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "FAM46B"] <- "TENT5B"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "FLJ34503"] <- "LINC02880"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "KIAA1024"] <- "MINAR1"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "MARC1"] <- "MTARC1"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "C9ORF3"] <- "AOPEP"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "LINC01279"] <- "CCDC80"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "LINC01468"] <- "LNCAROD"
  rownames(astrocytic_df)[rownames(astrocytic_df) == "DEPP1"] <- "C10orf10"
  
  
  astrocytic_df<-astrocytic_df[which(rownames(astrocytic_df) %in% names(TCGA_DATA)),]
  
  astrocytic_df$log2FC[astrocytic_df$log2FC < 0] <- "A1" 
  astrocytic_df$log2FC[astrocytic_df$log2FC != "A1"] <- "A2"
  RNAseq_genes<- as.vector(row.names(astrocytic_df))
  
  RNAseq_genes_astrocytic_df_list <- list(RNAseq_genes,astrocytic_df)
  
  return(RNAseq_genes_astrocytic_df_list)
}

# This function inserts A1/A2/progenitor/mature genes from Heiland paper into a vector.
# The function returns a list of: (1) vector of significant A1/A2/progenitor/mature genes (2) a df of all the genes with
# their A1/A2/progenitor/mature label.
Heiland_data_processing<- function(){
  Heiland_genes_df <- read.csv("1b+1e.csv")
  Heiland_genes_df <- data.frame(Heiland_genes_df[,-1], row.names = Heiland_genes_df[,1],check.names=FALSE)
  Heiland_genes_df<-Heiland_genes_df[ order(rownames(Heiland_genes_df)) , ,drop=F]
  Heiland_genes_df<-Heiland_genes_df[ which(row.names(Heiland_genes_df) %in% names(TCGA_DATA)), ,drop=F]
  names(Heiland_genes_df)[1]<- "label"
  Heiland_genes <- row.names(Heiland_genes_df)
  
  Heiland_genes_Heiland_df_list <- list(Heiland_genes,Heiland_genes_df)
  
  return(Heiland_genes_Heiland_df_list)
}

# This function inserts microglia associated genes from the microglia-dataset into a vector and
# returns the vector.
microglia_BRENDA_genes_processing<- function(){
  mocroglia_BRENDA_genes <- read.table("BRENDA_genes.txt", sep="\n", col.names=c("gene_name"))
  mocroglia_BRENDA_genes<-mocroglia_BRENDA_genes[ which(mocroglia_BRENDA_genes$gene_name %in% colnames(FPKM_DATA)[2:length((FPKM_DATA))]), ,drop=F]
  mocroglia_BRENDA_genes<-as.vector(mocroglia_BRENDA_genes$gene_name)
  return(mocroglia_BRENDA_genes)
}

# This function prepares data to a kaplan-meier analysis. 
# The data was downloaded in 2020. hence, the survival time is calculated by 2020-(year_of_diagnosis).
prep_for_KM <- function(km_data){
  km_data<-km_data[!(km_data$vital_status=="Not Reported"),]
  
  dead_patient_km<-km_data[km_data$vital_status=="Dead",]
  dead_patient_km<-dead_patient_km[!(dead_patient_km$days_to_death=="'--"),]
  dead_patient_km$days_to_death<- as.numeric(dead_patient_km$days_to_death)/365
  
  alive_patient_km<- km_data[km_data$vital_status=="Alive",]
  alive_patient_km$days_to_death<- 2020-(alive_patient_km$year_of_diagnosis)
  
  km_data<-rbind(dead_patient_km, alive_patient_km)
  km_data$year_of_diagnosis<-NULL
  km_data$days_to_death<- as.numeric(km_data$days_to_death)
  km_data$vital_status[km_data$vital_status == "Alive"] <- 0
  km_data$vital_status[km_data$vital_status == "Dead"] <- 1
  km_data$vital_status<- as.numeric(km_data$vital_status)
  
  km_data<- na.omit(km_data)  
  
  return(km_data)
}

# This function generates a kaplan-meier plot of a single given gene and checks if the plot is significant.
# the plot is not shown because we are interested in the p-val of the plot.
# the function returns the p-value of the plot.
run_km<- function(km_data, gene_name){
  if(gene_name %in% names(km_data)){
    km_data2<- dplyr::select(km_data, gene_name , days_to_death, vital_status)
    names(km_data2)[1]<-"gene"
    km_data2$gene<- as.numeric(km_data2$gene)
    km_data2$gene <- cut(km_data2$gene, breaks = c(-3,0,Inf),labels = c("low", "high"))
    km_data2<-km_data2[!(km_data2$gene=="mid"),]
    
    if(nrow(km_data2)==0){
      quit()
    }
    f1<-Surv(as.numeric(km_data2$days_to_death), as.numeric(km_data2$vital_status))
    fit = survfit(f1 ~ gene, data = km_data2)
    if(0) {
      ggsurvplot(fit, conf.int = TRUE,
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF", "#e978ff"),
                 pval= TRUE)
    }
    km_data2$vital_status<- as.numeric(km_data2$vital_status)
    fit2 <- surv_fit(Surv(days_to_death, vital_status) ~ gene, data = km_data2)
    pval<- surv_pvalue(fit2,method = "survdiff")[2]
    return(pval)
  }
  
}

# This funtion converts the survival column of  TCGA_DATA from numeric to categorical.
convert_num_surv_to_categorical_surv<- function(TCGA_DATA){
  
  TCGA_DATA$Survival..months.<- as.numeric(TCGA_DATA$Survival..months.)
  TCGA_DATA$new_survival<- NA
  TCGA_DATA<-subset(TCGA_DATA, TCGA_DATA$Survival..months.>0)
  
  for( i in 1:nrow(TCGA_DATA)){
    if(is.na(TCGA_DATA$Survival..months.[i])){
      TCGA_DATA$new_survival[i]<- "na"
    }else{
      if ((TCGA_DATA$Survival..months.[i] < 12) & (TCGA_DATA$Vital.status..1.dead.[i]=="1")){
        TCGA_DATA$new_survival[i]<- "less than 1"
        #TCGA_DATA$new_survival[i]<- TCGA_DATA$Survival..months.[i]
        
      }else{
        if (TCGA_DATA$Survival..months.[i] >= 12 & TCGA_DATA$Survival..months.[i] < 36 & TCGA_DATA$Vital.status..1.dead.[i]=="1"){
          TCGA_DATA$new_survival[i]<- "1-3"  
          #TCGA_DATA$new_survival[i]<- TCGA_DATA$Survival..months.[i]
          
        } else{
          if (TCGA_DATA$Survival..months.[i] >= 36){
            TCGA_DATA$new_survival[i]<- "more than 3"
            #TCGA_DATA$new_survival[i]<- TCGA_DATA$Survival..months.[i]
          } else{
            TCGA_DATA$new_survival[i]<- "na"
          } 
        }
      }
    }
  }
  TCGA_DATA<-cbind(TCGA_DATA[,which(colnames(TCGA_DATA) == "new_survival")], TCGA_DATA[,1:(ncol(TCGA_DATA)-1)])
  TCGA_DATA[,2]<- NULL
  names(TCGA_DATA)[1]<- "Survival..months."
  TCGA_DATA$Vital.status..1.dead.<- NULL
  
  return(TCGA_DATA)
}
