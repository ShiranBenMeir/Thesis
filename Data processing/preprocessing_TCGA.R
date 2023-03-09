## script i got from sahar and made a few changes. 
## The script reads the data that was downloaded from TCGA and creates "data_pre" data frame.
## data_pre- columns are the names the files (each name is the file_id + FPKM.txt and represent one patient).
##           rows are gene names in ENSG symbol (gene id).


setwd("C:/Users/shiran/Desktop/shiran/new CGA data/cart/zzz_all_unzip")
files <- list.files(path= getwd(),pattern = "\\.txt",
                    full.names = FALSE, recursive = FALSE,ignore.case = FALSE,include.dirs = FALSE)

data_pre<-NULL
data_pre[[1]] <- data.frame(read.table(files[1],sep="\t"))

for (i in 1:length(files)) {
  
  data_pre[[i]] <- data.frame(read.table(files[i],sep="\t"))
  colnames(data_pre[[i]]) <- c("GENEID",sub(".*/", "",files[i]))     # GENEID = ensembl ID; sub(".*/", "",files[i]) = file ID (taken from the manifest file)
  
  
}


if(1){
  temp<- data.frame()
  for (i in 1:(length(data_pre)-1)) {
    if(i==1){
      #all = F natural join
      temp<- Reduce(function(x,y) merge(x,y,by="GENEID",all=F) ,list(data_pre[[i]],data_pre[[i+1]]))
    }
    if(i>1){
      temp<- Reduce(function(x,y) merge(x,y,by="GENEID",all=F) ,list(temp,data_pre[[i+1]]))
    }
  }
  data_pre<-temp
  remove(temp)
  remove(i)
}

remove(files)
setwd("C:/Users/shiran/Documents/Rprojects/project2")