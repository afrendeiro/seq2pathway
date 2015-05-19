#function 1
addDescription <-
function(genome=c("hg19", "mm10", "mm9"),genevector){
if(missing(genome)){stop("Error: Please input genome type: 'mm10', 'mm9', or 'hg19'")}
if(genome == "mm10"){
   ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
   anno <- getBM(attributes=c("mgi_symbol", "mgi_description"), filters = "mgi_symbol", values = as.vector(genevector), mart = ensembl)
   anno <- anno[match(genevector, anno$mgi_symbol),]
   anno<-anno[is.na(anno$mgi_symbol)==FALSE,]
   }else if(genome == "mm9"){
   ensembl = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl',host="may2012.archive.ensembl.org")
   anno <- getBM(attributes=c("mgi_symbol", "mgi_description"), filters = "mgi_symbol", values = as.vector(genevector), mart = ensembl)
   anno <- anno[match(genevector, anno$mgi_symbol),]
   anno<-anno[is.na(anno$mgi_symbol)==FALSE,]
   }else if(genome == "hg19"){
   ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
   anno <- getBM(attributes=c("hgnc_symbol", "description"), filters = "hgnc_symbol", values = as.vector(genevector), mart = ensembl)
   anno <- anno[match(genevector, anno$hgnc_symbol),]
   anno<-anno[is.na(anno$hgnc_symbol)==FALSE,]
   }else{
   stop("Error: Please input genome type: 'mm10', 'mm9', or 'hg19'")
   }
   anno<-anno[is.na(anno[,1])==FALSE,]
   return(anno)
}

###function 2
cumulativerank <-
function(sampleExp, GeneID, GeneSet,logCheck,na.rm)
{

if(class(sampleExp)!="numeric"){sampleExp <-as.numeric(levels(sampleExp))[sampleExp]}
if(logCheck){if(max(sampleExp, na.rm=TRUE) > 20) {sampleExp <- log2(sampleExp)}}
if(length(sampleExp)!=length(GeneID)){stop("Error: GeneID information is missing or not correct!")}

N <- length(GeneID)
GeneID_NaInSet <- GeneID[which(!GeneID %in% GeneSet)]

# Step 1
    rankedExp <- rank(sampleExp, na.last="keep")
    rankscore <- rankedExp

# Step 2:

x1 <- which(GeneID %in% GeneSet)
x2 <- which(GeneID %in% GeneID_NaInSet)

ST <- sum(rankscore[x1], na.rm=TRUE)/length(x1)
SN <- sum(rankscore[x2], na.rm=TRUE)/length(x2)
y <- sum(ST, -SN, na.rm=na.rm)


return(y)
}

###function3
cumulativerank_EmpiricalP <-
function(sampleExp, GeneID, GeneSet, logCheck,cumulativerank,B,na.rm)
{
data(gencode_coding,package="seq2pathway.data")
if(class(sampleExp)!="numeric"){sampleExp <-as.numeric(levels(sampleExp))[sampleExp]}
if(logCheck){if(max(sampleExp, na.rm=TRUE) > 20) {sampleExp <- log2(sampleExp)}}
if(length(sampleExp)!=length(GeneID)){stop("Error: GeneID information is missing or not correct!")}

N <- length(GeneID)
count<-0
if(is.na(cumulativerank)==FALSE){

#Step 1: Calculation of weighted rank of gene expression
    rankedExp <- rank(sampleExp, na.last="keep")
    rankscore <- rankedExp

for(test in 1:B){
GeneID<-sample(gencode_coding,N)
GeneID_NaInSet <- GeneID[which(!GeneID %in% GeneSet)]

# Step 2:
x1 <- which(GeneID %in% GeneSet)
x2 <- which(GeneID %in% GeneID_NaInSet)


ST <- sum(rankscore[x1], na.rm=TRUE)/length(x1)
SN <- sum(rankscore[x2], na.rm=TRUE)/length(x2)
y <- sum(ST, -SN, na.rm=na.rm)

if(y>=cumulativerank&is.na(y)==FALSE){count=count+1}
}#B's sampling loop ends

Pvalue<-count/B
}#end of if cumulativerank is not equal to NA
else{Pvalue=NA}

return(Pvalue)
}

#Function 4
FAIME <-
function(sampleExp, GeneID, GeneSet, alpha, logCheck,na.rm)
{
if(class(sampleExp)!="numeric"){sampleExp <-as.numeric(levels(sampleExp))[sampleExp]}
if(logCheck){if(max(sampleExp, na.rm=TRUE) > 20) {sampleExp <- log2(sampleExp)}}
if(length(sampleExp)!=length(GeneID)){stop("Error: GeneID information is missing or not correct!")}

N <- length(GeneID)
GeneID_NaInSet <- GeneID[which(!GeneID %in% GeneSet)]

#Step 1: Calculation of weighted rank of gene expression
    rankedExp <- rank(sampleExp, na.last="keep")
    rankscore <- rankedExp*exp((rankedExp/N*alpha)-alpha)

# Step 2:

x1 <- which(GeneID %in% GeneSet)
x2 <- which(GeneID %in% GeneID_NaInSet)

ST <- sum(rankscore[x1], na.rm=TRUE)/length(x1)
SN <- sum(rankscore[x2], na.rm=TRUE)/length(x2)
y <- sum(ST, -SN, na.rm=na.rm)

return(y)
}


#function 5
FAIME_EmpiricalP <-
function(sampleExp, GeneID, GeneSet, alpha, logCheck,FAIME,B,na.rm)
{

data(gencode_coding,package="seq2pathway.data")
if(class(sampleExp)!="numeric"){sampleExp <-as.numeric(levels(sampleExp))[sampleExp]}
if(logCheck){if(max(sampleExp, na.rm=TRUE) > 20) {sampleExp <- log2(sampleExp)}}
if(length(sampleExp)!=length(GeneID)){stop("Error: GeneID information is missing or not correct!")}

N <- length(GeneID)
count<-0
if(is.na(FAIME)==FALSE){

#Step 1: Calculation of weighted rank of gene expression
    rankedExp <- rank(sampleExp, na.last="keep")
    rankscore <- rankedExp*exp((rankedExp/N*alpha)-alpha)

for(test in 1:B){
GeneID<-sample(as.vector(gencode_coding),N)
GeneID_NaInSet <- GeneID[which(!GeneID %in% GeneSet)]

# Step 2:
x1 <- which(GeneID %in% GeneSet)
x2 <- which(GeneID %in% GeneID_NaInSet)


ST <- sum(rankscore[x1], na.rm=TRUE)/length(x1)
SN <- sum(rankscore[x2], na.rm=TRUE)/length(x2)
y <- sum(ST, -SN, na.rm=na.rm)

if(y>=FAIME&is.na(y)==FALSE){count=count+1}
}#B's sampling loop ends

Pvalue<-count/B
}#end of if FAIME is not equal to NA
else{Pvalue=NA}

return(Pvalue)
}


#Function 6
FisherTest_GO_BP_MF_CC <-
function(gs,genome=c("hg38","hg19","mm10","mm9"),min_Intersect_Count=5){
if(missing(min_Intersect_Count)){min_Intersect_Count=5}

####load GP pathway information
data(GO_BP_list,package="seq2pathway.data")
data(GO_MF_list,package="seq2pathway.data")
data(GO_CC_list,package="seq2pathway.data")

data(Des_BP_list,package="seq2pathway.data")
data(Des_MF_list,package="seq2pathway.data")
data(Des_CC_list,package="seq2pathway.data")

#######check genome
if(missing(genome)){genome="hg19"}
#####load  GO_GENCODE_databse_intersect_gene information
if(genome=="hg38"){
data(GO_GENCODE_df_hg_v20,package="seq2pathway.data")
GO_GENCODE_df<-GO_GENCODE_df_hg_v20
}else if(genome=="hg19"){
data(GO_GENCODE_df_hg_v19,package="seq2pathway.data")
GO_GENCODE_df<-GO_GENCODE_df_hg_v19
}else if(genome=="mm10"){
data(GO_GENCODE_df_mm_vM3,package="seq2pathway.data")
GO_GENCODE_df<-GO_GENCODE_df_mm_vM3
}else if(genome=="mm9"){
data(GO_GENCODE_df_mm_vM1,package="seq2pathway.data")
GO_GENCODE_df<-GO_GENCODE_df_mm_vM1
}

GO_FisherTest_result<-list()

####################GO:BP
###background counts
Para_GenesetsGene<-length(intersect(unique(toupper(unlist(GO_BP_list))), union(unique(toupper(GO_GENCODE_df$gene_name)),toupper(unlist(gs)))))
###gene vector annotated from GENCODE intersect with GO_BP_list
Para_PeakGene<-length(intersect(unique(toupper(unlist(GO_BP_list))),unique(toupper(unlist(gs)))))


mdat <- matrix(, nrow = length(GO_BP_list), ncol = 9, byrow = TRUE,
               dimnames = list(c(1:length(GO_BP_list)),
               c("GOID","Description","Fisher_Pvalue", "Fisher_odds","FDR","Intersect_Count","GO_gene_inBackground","GO_gene_raw_Count","Intersect_gene")))
mdat<-as.data.frame(mdat)
for(i in 1:length(GO_BP_list)){

pathwaygene<-length(intersect(toupper(GO_BP_list[[i]]),union(unique(toupper(GO_GENCODE_df$gene_name)),toupper(unlist(gs)))))

a<-length(intersect(toupper(GO_BP_list[[i]]),unique(toupper(unlist(gs)))))
b<-Para_PeakGene-a
c<-pathwaygene-a
d<-Para_GenesetsGene-a-b-c

Con<-matrix(c(a, b, c, d),nrow = 2, byrow =TRUE)
mdat[i,1]<-names(GO_BP_list)[i]
mdat[i,2]<-as.character(Des_BP_list[which(names(Des_BP_list)==names(GO_BP_list)[i])])
mdat[i,3]<-fisher.test(Con)$p.value
mdat[i,4]<-fisher.test(Con)$estimate
mdat[i,6]<-a
mdat[i,7]<-pathwaygene
mdat[i,8]<-length(GO_BP_list[[i]])
mdat[i,9]<-paste(intersect(toupper(GO_BP_list[[i]]),unique(toupper(unlist(gs)))),collapse=" ")
}
mdat<-mdat[mdat$Intersect_Count>=min_Intersect_Count,]
mdat$FDR <- p.adjust(mdat[,c("Fisher_Pvalue")],method="BH")
mdat<-mdat[order(mdat$FDR),]
row.names(mdat)<-NULL
GO_FisherTest_result[[1]]<-mdat
rm(mdat)
names(GO_FisherTest_result)[1]<-c("GO_BP")

####################GO:CC
###background counts
Para_GenesetsGene<-length(intersect(unique(toupper(unlist(GO_CC_list))), union(unique(toupper(GO_GENCODE_df$gene_name)),toupper(unlist(gs)))))
###gene vector annotated from GENCODE intersect with GO_CC_list
Para_PeakGene<-length(intersect(unique(toupper(unlist(GO_CC_list))),unique(toupper(unlist(gs)))))


mdat <- matrix(, nrow = length(GO_CC_list), ncol = 9, byrow = TRUE,
               dimnames = list(c(1:length(GO_CC_list)),
                               c("GOID","Description","Fisher_Pvalue", "Fisher_odds","FDR","Intersect_Count","GO_gene_inBackground","GO_gene_raw_Count","Intersect_gene")))
mdat<-as.data.frame(mdat)
for(i in 1:length(GO_CC_list)){

pathwaygene<-length(intersect(toupper(GO_CC_list[[i]]),union(unique(toupper(GO_GENCODE_df$gene_name)),toupper(unlist(gs)))))

a<-length(intersect(toupper(GO_CC_list[[i]]),unique(toupper(unlist(gs)))))
b<-Para_PeakGene-a
c<-pathwaygene-a
d<-Para_GenesetsGene-a-b-c

Con<-matrix(c(a, b, c, d),nrow = 2, byrow =TRUE)
mdat[i,1]<-names(GO_CC_list)[i]
mdat[i,2]<-Des_CC_list[which(names(Des_CC_list)==names(GO_CC_list)[i])]
mdat[i,3]<-fisher.test(Con)$p.value
mdat[i,4]<-fisher.test(Con)$estimate
mdat[i,6]<-a
mdat[i,7]<-pathwaygene
mdat[i,8]<-length(GO_CC_list[[i]])
mdat[i,9]<-paste(intersect(toupper(GO_CC_list[[i]]),unique(toupper(unlist(gs)))),collapse=" ")
}
mdat<-mdat[mdat$Intersect_Count>=min_Intersect_Count,]
mdat$FDR <- p.adjust(mdat[,c("Fisher_Pvalue")],method="BH")
mdat<-mdat[order(mdat$FDR),]
row.names(mdat)<-NULL
GO_FisherTest_result[[2]]<-mdat
rm(mdat)
names(GO_FisherTest_result)[2]<-c("GO_CC")

####################GO:MF
###background counts
Para_GenesetsGene<-length(intersect(unique(toupper(unlist(GO_MF_list))),union(unique(toupper(GO_GENCODE_df$gene_name)),toupper(unlist(gs)))))
###gene vector annotated from GENCODE intersect with GO_MF_list
Para_PeakGene<-length(intersect(unique(toupper(unlist(GO_MF_list))),unique(toupper(unlist(gs)))))


mdat <- matrix(, nrow = length(GO_MF_list), ncol = 9, byrow = TRUE,
               dimnames = list(c(1:length(GO_MF_list)),
                               c("GOID","Description","Fisher_Pvalue", "Fisher_odds","FDR","Intersect_Count","GO_gene_inBackground","GO_gene_raw_Count","Intersect_gene")))
mdat<-as.data.frame(mdat)
for(i in 1:length(GO_MF_list)){

pathwaygene<-length(intersect(toupper(GO_MF_list[[i]]),union(unique(toupper(GO_GENCODE_df$gene_name)),toupper(unlist(gs)))))

a<-length(intersect(toupper(GO_MF_list[[i]]),unique(toupper(unlist(gs)))))
b<-Para_PeakGene-a
c<-pathwaygene-a
d<-Para_GenesetsGene-a-b-c

Con<-matrix(c(a, b, c, d),nrow = 2, byrow =TRUE)
mdat[i,1]<-names(GO_MF_list)[i]
mdat[i,2]<-Des_MF_list[which(names(Des_MF_list)==names(GO_MF_list)[i])]
mdat[i,3]<-fisher.test(Con)$p.value
mdat[i,4]<-fisher.test(Con)$estimate
mdat[i,6]<-a
mdat[i,7]<-pathwaygene
mdat[i,8]<-length(GO_MF_list[[i]])
mdat[i,9]<-paste(intersect(toupper(GO_MF_list[[i]]),unique(toupper(unlist(gs)))),collapse=" ")
}
mdat<-mdat[mdat$Intersect_Count>=min_Intersect_Count,]
mdat$FDR <- p.adjust(mdat[,c("Fisher_Pvalue")],method="BH")
mdat<-mdat[order(mdat$FDR),]
row.names(mdat)<-NULL
GO_FisherTest_result[[3]]<-mdat
rm(mdat)
names(GO_FisherTest_result)[3]<-c("GO_MF")

print("Fisher's exact test done")
return(GO_FisherTest_result)
}

#function 7
FisherTest_MsigDB <-
function(gsmap,gs,genome=c("hg38","hg19","mm10","mm9"),min_Intersect_Count=5){
if(missing(min_Intersect_Count)){min_Intersect_Count=5}

if(class(gsmap)!="GSA.genesets"){stop("Error: gsmap should be GAS.genesets object")}
#######check genome
if(missing(genome)){genome="hg19"}
if(genome=="hg38"){
data(Msig_GENCODE_df_hg_v20,package="seq2pathway.data")
Msig_GENCODE_df<-Msig_GENCODE_df_hg_v20
}else if(genome=="hg19"){
data(Msig_GENCODE_df_hg_v19,package="seq2pathway.data")
Msig_GENCODE_df<-Msig_GENCODE_df_hg_v19
}else if(genome=="mm10"){
data(Msig_GENCODE_df_mm_vM3,package="seq2pathway.data")
Msig_GENCODE_df<-Msig_GENCODE_df_mm_vM3
}else if(genome=="mm9"){
data(Msig_GENCODE_df_mm_vM1,package="seq2pathway.data")
Msig_GENCODE_df<-Msig_GENCODE_df_mm_vM1
}


###background counts
Para_GenesetsGene<-length(intersect(unique(toupper(unlist(gsmap$genesets))), union(unique(toupper(Msig_GENCODE_df$gene_name)),unique(toupper(unlist(gs))))))

###gene vector in Msigdb
Para_PeakGene<-length(intersect(unique(toupper(unlist(gsmap$genesets))),unique(toupper(unlist(gs)))))



mdat <- matrix(, nrow = length(gsmap$geneset.names), ncol = 9, byrow = TRUE,
               dimnames = list(c(1:length(gsmap$geneset.names)),
               c("GeneSet","Description","Fisher_Pvalue","Fisher_odds","FDR","Intersect_Count","MsigDB_gene_inBackground","MsigDB_gene_raw_Count","Intersect_gene")))
mdat<-as.data.frame(mdat)
for(i in 1:length(gsmap$geneset.names)){

pathwaygene<-length(intersect(toupper(gsmap$genesets[[i]]),union(unique(toupper(Msig_GENCODE_df$gene_name)),toupper(unlist(gs)))))

a<-length(intersect(toupper(gsmap$genesets[[i]]),unique(toupper(unlist(gs)))))
b<-Para_PeakGene-a
c<-pathwaygene-a
d<-Para_GenesetsGene-a-b-c

Con<-matrix(c(a, b, c, d),nrow = 2, byrow =TRUE)
mdat[i,1]<-gsmap$geneset.names[i]
mdat[i,2]<-gsmap$geneset.descriptions[i]
mdat[i,3]<-fisher.test(Con)$p.value
mdat[i,4]<-fisher.test(Con)$estimate
mdat[i,6]<-a
mdat[i,7]<-pathwaygene
mdat[i,8]<-length(gsmap$genesets[[i]])
mdat[i,9]<-paste(intersect(toupper(gsmap$genesets[[i]]),unique(toupper(unlist(gs)))),collapse=" ")
}

mdat<-mdat[mdat$Intersect_Count>=min_Intersect_Count,]
mdat$FDR <- p.adjust(mdat[,c("Fisher_Pvalue")],method="BH")
mdat<-mdat[order(mdat$FDR),]
row.names(mdat)<-NULL

print("Fisher's exact test done")
return(mdat)
}

#function 8
KSrank <-
function(sampleExp, GeneID, GeneSet, alpha, logCheck)
{
if(class(sampleExp)!="numeric"){sampleExp <-as.numeric(levels(sampleExp))[sampleExp]}
if(logCheck){if(max(sampleExp, na.rm=TRUE) > 20) {sampleExp <- log2(sampleExp)}}
if(length(sampleExp)!=length(GeneID)){stop("Error: GeneID information is missing or not correct!")}

N <- length(GeneID)
GeneID_NaInSet <- GeneID[which(!GeneID %in% GeneSet)]

# Step 1 (Supporting Figure 2 in the Text S1): Calculation of weighted rank of gene expression -----
## Ranked by score,  the lowest to highest. Therefore, the up-regulated genes get the higher weighted score
#Only dir="high"
    rankedExp <- rank(sampleExp, na.last="keep")
    rankscore <- rankedExp*exp((rankedExp/N*alpha)-alpha)

# Step 2:

x1 <- which(GeneID %in% GeneSet)
x2 <- which(GeneID %in% GeneID_NaInSet)

## WFA algorithm, KS test is emploied
   if(length(x1)>0 & length(x2)>0) {
         tmp <- ks.test(rankscore[x1], rankscore[x2], alternative="two.sided")
         y <- tmp$statistic

   }else{y=NA}

return(y)
}


#function 9
KSrank_EmpiricalP <-
function(sampleExp, GeneID, GeneSet, alpha, logCheck,KSrank,B)
{
data(gencode_coding,package="seq2pathway.data")
if(class(sampleExp)!="numeric"){sampleExp <-as.numeric(levels(sampleExp))[sampleExp]}
if(logCheck){if(max(sampleExp, na.rm=TRUE) > 20) {sampleExp <- log2(sampleExp)}}
if(length(sampleExp)!=length(GeneID)){stop("Error: GeneID information is missing or not correct!")}

N <- length(GeneID)
count<-0
if(is.na(KSrank)==FALSE){

#Step 1: Calculation of weighted rank of gene expression
    rankedExp <- rank(sampleExp, na.last="keep")
    rankscore <- rankedExp*exp((rankedExp/N*alpha)-alpha)

for(test in 1:B){
GeneID<-sample(gencode_coding,N)
GeneID_NaInSet <- GeneID[which(!GeneID %in% GeneSet)]

# Step 2:
x1 <- which(GeneID %in% GeneSet)
x2 <- which(GeneID %in% GeneID_NaInSet)


## WFA algorithm, KS test is emploied
   if(length(x1)>0 & length(x2)>0) {
         tmp <- ks.test(rankscore[x1], rankscore[x2], alternative="two.sided")
         y <- tmp$statistic

   }else{y=NA}

if(y>=KSrank&is.na(y)==FALSE){count=count+1}
}#B's sampling loop ends

Pvalue<-count/B
}#end of if KSrank is not equal to NA
else{Pvalue=NA}

return(Pvalue)
}


#function 10
Normalize_F <-
function(input){
Y.high1 <- apply(input, 2, scale)
rownames(Y.high1)<-rownames(input)
colnames(Y.high1)<-c(paste(colnames(input),"_Normalized",sep=""))
dim(Y.high1)
head(Y.high1)
cp2_FAIME<-as.data.frame(Y.high1)
print("Normalization.........done")
return(cp2_FAIME)
}

#function 11
Peak_Gene_Collapse <-
function(input,collapsemethod){
options(warn=-1)

#gene symbol or entrezID
GeneGroup<-as.vector(toupper(input[,2]))

#rowID must be distint,rename rownames of input
rownames(input)<-1:nrow(input)
rowID<-rownames(input)

#just for single sample file, otherwise, collpaseRows function could not work
input$index=0

#drop the first two columns, peak id and gene information
dataforCollapse<-input[,3:dim(input)[2]]

#main function
Collapse_Gene<-collapseRows(dataforCollapse, GeneGroup, rowID,  method=collapsemethod)

#main output of collapsed matrix, the rownames are gene symbol/entrezID
PeakGene_MaxMean<-Collapse_Gene$datETcollapsed

##get the reserved rownumber
rw<-as.numeric(Collapse_Gene$group2row[,2])

#remove the dummy index column, and the gene symbol column
selected_Peak_Info<-input[rw,c(1,3:(dim(input)[2]-1))]#remove the dummy index column

##gene symbol or entrezID as rownames
rownames(selected_Peak_Info)<-rownames(PeakGene_MaxMean)

###dummmy yummy
selected_Peak_Info1<-data.frame(selected_Peak_Info)
rownames(selected_Peak_Info1)<-rownames(selected_Peak_Info)
names(selected_Peak_Info1)<-names(selected_Peak_Info)
print("Peak_Gene_Collapse....... done")
return(selected_Peak_Info1)
}

#function 12
rungene2pathway <-
function(dat, gsmap, alpha=5, logCheck=FALSE, method=c("FAIME","KS-rank","cumulative-rank"),na.rm=FALSE)
{
if(missing(method)){method="FAIME"}
if(missing(alpha)){alpha=5}
if(missing(na.rm)){na.rm=FALSE}
if(missing(logCheck)){logCheck=FALSE}

if(! method %in% c("FAIME", "KS-rank", "cumulative-rank")) {stop("Error: method should be a value of 'FAIME', 'KS-rank', or 'cumulative-rank'!")}
if(is.data.frame(dat)==FALSE & is.matrix(dat)==FALSE){stop("Error: input should be a data frame or a matrix")}

if(class(gsmap)=="GSA.genesets"){
   seeds <- gsmap$geneset.names
   res <- matrix(nrow=length(seeds), ncol=ncol(dat))
   for(i in 1:length(seeds)){
           for(j in 1:ncol(dat)){
              if(method=="FAIME"){res[i,j] <- FAIME(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), alpha=alpha, logCheck=logCheck,na.rm=na.rm)
              }else if(method=="KS-rank"){res[i,j] <- KSrank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), alpha=alpha, logCheck=logCheck)
              }else if(method=="cumulative-rank"){res[i,j] <- cumulativerank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), logCheck=logCheck,na.rm=na.rm)
              }

           }#j loop
       }#i loop
   }else if(class(gsmap)=="list"){
    if(is.null(names(gsmap))){stop ("please give the names of gsmap as a list")}
    seeds <- names(gsmap)
    res <- matrix(nrow=length(seeds), ncol=ncol(dat))
    for(i in 1:length(seeds)){
           for(j in 1:ncol(dat)){
               if(method=="FAIME"){res[i,j] <- FAIME(sampleExp=dat[,j],GeneID=toupper(rownames(dat)),GeneSet= toupper(gsmap[[i]]), alpha=alpha, logCheck=logCheck,na.rm=na.rm)
               }else if(method=="KS-rank"){res[i,j] <- KSrank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap[[i]]), alpha=alpha, logCheck=logCheck)
               }else if(method=="cumulative-rank"){res[i,j] <- cumulativerank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap[[i]]), logCheck=logCheck,na.rm=na.rm)
               }

          }#j loop
       }#i loop
    }#class loop

rownames(res) <- seeds
colnames(res) <- c(paste(colnames(dat),"2pathscore",sep=""))
print("gene2pathway calculates score....... done")
return(res)
}


#function 13
rungene2pathway_EmpiricalP <-
function(dat, gsmap, alpha=5, logCheck=FALSE, method=c("FAIME","KS-rank","cumulative-rank"),B=100,na.rm=FALSE)
{
if(missing(method)){method="FAIME"}
if(missing(alpha)){alpha=5}
if(missing(na.rm)){na.rm=FALSE}
if(missing(B)){B=100}
if(missing(logCheck)){logCheck=FALSE}

if(! method %in% c("FAIME", "KS-rank", "cumulative-rank")) {stop("Error: method should be a value of 'FAIME', 'KS-rank', or 'cumulative-rank'!")}
if(is.data.frame(dat)==FALSE & is.matrix(dat)==FALSE){stop("Error: input should be a data frame or a matrix")}

if(class(gsmap)=="GSA.genesets"){
   seeds <- gsmap$geneset.names
   res <- matrix(nrow=length(seeds), ncol=ncol(dat))
   for(i in 1:length(seeds)){
           for(j in 1:ncol(dat)){
              if(method=="FAIME"){res[i,j] <- FAIME(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), alpha=alpha,logCheck=logCheck,na.rm=na.rm)
              }else if(method=="KS-rank"){res[i,j] <- KSrank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), alpha=alpha, logCheck=logCheck)
              }else if(method=="cumulative-rank"){res[i,j] <- cumulativerank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), logCheck=logCheck,na.rm=na.rm)
              }

           }#j loop
       }#i loop
   }else if(class(gsmap)=="list"){
    if(is.null(names(gsmap))){stop ("please give the names of gsmap as a list")}
    seeds <- names(gsmap)
    res <- matrix(nrow=length(seeds), ncol=ncol(dat))
    for(i in 1:length(seeds)){
           for(j in 1:ncol(dat)){
               if(method=="FAIME"){res[i,j] <- FAIME(sampleExp=dat[,j],GeneID=toupper(rownames(dat)),GeneSet= toupper(gsmap[[i]]), alpha=alpha, logCheck=logCheck,na.rm=na.rm)
               }else if(method=="KS-rank"){res[i,j] <- KSrank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap[[i]]), alpha=alpha, logCheck=logCheck)
               }else if(method=="cumulative-rank"){res[i,j] <- cumulativerank(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap[[i]]), logCheck=logCheck,na.rm=na.rm)
               }

          }#j loop
       }#i loop
    }#class loop

rownames(res) <- seeds
colnames(res) <- c(paste(colnames(dat),"2pathscore",sep=""))

#######Just progress bar
pb <- txtProgressBar(min = 0, max = length(seeds), style = 3)
#######Just progress bar


if(class(gsmap)=="GSA.genesets"){
   seeds <- gsmap$geneset.names
   res_p <- matrix(nrow=length(seeds), ncol=ncol(dat))
   for(i in 1:length(seeds)){
           #if(i%%100 ==0){print(paste("still working..........",Sys.time())) }
           setTxtProgressBar(pb, i)
           for(j in 1:ncol(dat)){
              if(method=="FAIME"){res_p[i,j] <- FAIME_EmpiricalP(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), alpha=alpha, logCheck=logCheck,FAIME=res[i,j],B=B,na.rm=na.rm)
              }else if(method=="KS-rank"){res_p[i,j] <- KSrank_EmpiricalP(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), alpha=alpha, logCheck=logCheck,KSrank=res[i,j],B=B)
              }else if(method=="cumulative-rank"){res_p[i,j] <-cumulativerank_EmpiricalP(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap$genesets[[i]]), logCheck=logCheck,cumulativerank=res[i,j],B=B,na.rm=na.rm)
              }
           }#j loop
       }#i loop

   }else if(class(gsmap)=="list"){
    if(is.null(names(gsmap))){stop ("please give the names of gsmap as a list")}
    seeds <- names(gsmap)
    res_p <- matrix(nrow=length(seeds), ncol=ncol(dat))
    for(i in 1:length(seeds)){
           #if(i%%100 ==0){print(paste("still working..........",Sys.time()))}
	   setTxtProgressBar(pb, i)
           for(j in 1:ncol(dat)){
               if(method=="FAIME"){res_p[i,j]<- FAIME_EmpiricalP(sampleExp=dat[,j],GeneID=toupper(rownames(dat)),GeneSet= toupper(gsmap[[i]]), alpha=alpha, logCheck=logCheck,FAIME=res[i,j],B=B,na.rm=na.rm)
               }else if(method=="KS-rank"){res_p[i,j] <- KSrank_EmpiricalP(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap[[i]]), alpha=alpha, logCheck=logCheck,KSrank=res[i,j],B=B)
               }else if(method=="cumulative-rank"){res_p[i,j] <-cumulativerank_EmpiricalP(sampleExp=dat[,j], GeneID=toupper(rownames(dat)), GeneSet=toupper(gsmap[[i]]), logCheck=logCheck,cumulativerank=res[i,j],B=B,na.rm=na.rm)
               }

          }#j loop
       }#i loop
    }#class loop


rownames(res_p) <- seeds
colnames(res_p) <- c(paste(colnames(dat),"2pathscore_Pvalue",sep=""))
#######Just progress bar
close(pb)
print("pathwayscore Empirical Pvalue calculation..........done")
return(res_p)

}

#function 14
runseq2gene <-
    function(inputfile,search_radius=150000,promoter_radius=200,promoter_radius2=100,genome=c("hg38","hg19","mm10","mm9"),
           adjacent=FALSE,SNP= FALSE,PromoterStop=FALSE,NearestTwoDirection=TRUE,UTR3=FALSE){
    
    options(scipen=999)
    if(missing(inputfile)){stop("please give the input file")}
    if(! class(inputfile) %in% c("data.frame","GRanges")){stop("please check the format of input file")}
    
    ###default parameters if missing
    if(missing(genome)){genome= "hg19"}
    if(missing(search_radius)){search_radius = 150000}
    if(missing(promoter_radius)){promoter_radius = 200}
    if(missing(promoter_radius2)){promoter_radius2 = 100}
    if(missing(SNP)){SNP = "False"}
    if(missing(adjacent)){adjacent= "False"}
    if(missing(PromoterStop)){PromoterStop= "False"}
    if(missing(NearestTwoDirection)){NearestTwoDirection= "True"}
    if(missing(UTR3)){UTR3= "False"}
    
    if(SNP %in% c("T","TRUE","True",TRUE)){SNP= "True"}
    if(SNP %in% c("F","FALSE","False",FALSE)){SNP= "False"}
    
    if(PromoterStop %in% c("T","TRUE","True",TRUE)){PromoterStop= "True"}
    if(PromoterStop %in% c("F","FALSE","False",FALSE)){PromoterStop= "False"}
    
    if(NearestTwoDirection %in% c("T","TRUE","True",TRUE)){NearestTwoDirection= "True"}
    if(NearestTwoDirection %in% c("F","FALSE","False",FALSE)){NearestTwoDirection= "False"}
    
    if(UTR3 %in% c("T","TRUE","True",TRUE)){UTR3= "True"}
    if(UTR3 %in% c("F","FALSE","False",FALSE)){UTR3= "False"}
    
    if(adjacent %in% c("T","TRUE","True",TRUE)){adjacent= "True"}
    if(adjacent %in% c("F","FALSE","False",FALSE)){adjacent= "False"}
    if(adjacent== "True"){search_radius =0}
    
    if(length(genome>1)){genome=genome[1]}
    
    ### assign the path of main function
    path<-paste(system.file(package="seq2pathway"),"/scripts/Function_PeakMutationAnnotation_GENCODE_05182015.py",sep="/")
    
    ###wrap invoke file
    name<-paste("inputfile",gsub(":","_",gsub("-","",gsub(" ","_",Sys.time()))),"seq2gene_log.py",sep="_")
    
    tmpinfile = tempfile()
    if(class(inputfile)=="data.frame"){
      if(ncol(inputfile)<4){stop("please check the format of the input data, some column is missing")}
      write.table(inputfile, file=tmpinfile,sep="\t",quote=FALSE,row.names = FALSE)}
    if(class(inputfile)=="GRanges"){
      test<-as.data.frame(inputfile)
      if(ncol(test)<6){stop("please check the format of the input data, some column is missing")}
      if(ncol(test)==6){write.table(test[,c(6,1,2,3)], file=tmpinfile,sep="\t",quote=FALSE,row.names = FALSE)}
      if(ncol(test)>6){write.table(test[,c(6,1,2,3,7)], file=tmpinfile,sep="\t",quote=FALSE,row.names = FALSE)}
    }
    tmpinfile = gsub("\\","/",tmpinfile,fixed =TRUE)
    tmpoutfile = tempfile()
    write.table(NULL, file=tmpoutfile,sep="\t",quote=FALSE,row.names = FALSE)
    tmpoutfile = gsub("\\","/",tmpoutfile,fixed =TRUE)

    tmpoutfile_UTR3 = tempfile()
    write.table(NULL, file=tmpoutfile_UTR3,sep="\t",quote=FALSE,row.names = FALSE)
    tmpoutfile_UTR3 = gsub("\\","/",tmpoutfile_UTR3,fixed =TRUE)
    
    tmp_ref_file = tempfile()
    write.table(NULL, file=tmp_ref_file,sep="\t",quote=FALSE,row.names = FALSE)
    tmp_ref_file = gsub("\\","/",tmp_ref_file,fixed =TRUE)
    
    sink(paste(tempdir(),"\\",name,sep=""))
    ###fixed headers import modules
    cat("import sys, string, math, shutil, math, os, gzip, time, glob, multiprocessing",sep="\n")
    cat("from shutil import rmtree",sep="\n")
    cat("from datetime import datetime",sep="\n")
    cat("from bisect import *",sep="\n")
    cat("",sep="\n")
    ###import our function module
    cat("import imp",sep="\n")
    cat("imp.load_source('Function_PeakMutationAnnotation_GENCODE_05182015',")
    cat("'", path, "')",sep="")
    cat("",sep="\n")
    cat("from Function_PeakMutationAnnotation_GENCODE_05182015 import FindPeakMutation",sep="\n")
    cat("",sep="\n")
    
    ####write parameters
    #cat(paste("inputpath=","'",inputpath,"/'",sep=""),sep="\n")
    cat(paste("inputfile=","'",tmpinfile,"'",sep=""),sep="\n")
    #cat(paste("outputpath=","'",outputpath,"/'",sep=""),sep="\n")
    cat(paste("outputfile=","'",tmpoutfile,"'",sep=""),sep="\n")
    cat(paste("outputfileUTR3=","'",tmpoutfile_UTR3,"'",sep=""),sep="\n")
    cat(paste("tmp_ref_file=","'",tmp_ref_file,"'",sep=""),sep="\n")
    cat(paste("search_radius=",search_radius,sep=""),sep="\n")
    cat(paste("promoter_radius=",promoter_radius,sep=""),sep="\n")
    cat(paste("promoter_radius2=",promoter_radius2,sep=""),sep="\n")  
    cat(paste("genome=","'",genome,"'",sep=""),sep="\n")
    cat(paste("adjacent=",adjacent,sep=""),sep="\n")
    cat(paste("pwd=","'",system.file(package="seq2pathway.data"),"/extdata/'",sep=""),sep="\n")
    cat(paste("SNP=",SNP,sep=""),sep="\n")
    cat(paste("PromoterStop=",PromoterStop,sep=""),sep="\n") 
    cat(paste("NearestTwoDirection=",NearestTwoDirection,sep=""),sep="\n")
    cat(paste("UTR3=",UTR3,sep=""),sep="\n")
    cat("FindPeakMutation(inputfile,outputfile,outputfileUTR3,tmp_ref_file,search_radius,promoter_radius,promoter_radius2,genome,adjacent,pwd,SNP,PromoterStop,NearestTwoDirection,UTR3)")
    sink()
    
    #invoke python
    command <- paste("C:/Python27/python ", tempdir(),"\\",name,sep="")
    response <- system(command, intern=TRUE)
    anno_result<-read.table(file=tmpoutfile_UTR3,header=TRUE,sep="\t")
    seq2gene_result=list()
    seq2gene_result[[1]]<-anno_result
    names(seq2gene_result)[1]<-"seq2gene_FullResult"
    seq2gene_result[[2]]<-anno_result[anno_result$source=="protein_coding",]
    names(seq2gene_result)[2]<-"seq2gene_CodingGeneOnlyResult"  
    print(response)
    return(seq2gene_result)
  }



#function 15
runseq2pathway<-function(inputfile,
         search_radius=150000,promoter_radius=200,promoter_radius2=100,
         genome=c("hg38","hg19","mm10","mm9"),adjacent=FALSE,SNP= FALSE,PromoterStop=FALSE,NearestTwoDirection=TRUE,UTR3=FALSE,
         DataBase=c("GOterm"),FAIMETest=FALSE,FisherTest=TRUE,
         collapsemethod=c("MaxMean","function","ME","maxRowVariance","MinMean","absMinMean","absMaxMean","Average"),
         alpha=5,logCheck=FALSE,B=100,na.rm=FALSE,min_Intersect_Count=5){
options(warn=-1)

if(missing(inputfile)){stop("please give the input file")}
if(! class(inputfile) %in% c("data.frame","GRanges")){stop("please check the format of input file")}

###default parameters if missing
if(missing(DataBase)){DataBase="GOterm"}
if(missing(FAIMETest)){FAIMETest=FALSE}
if(missing(FisherTest)){FisherTest=TRUE}
if(missing(genome)){genome= "hg19"}
if(missing(search_radius)){search_radius = 150000}
if(missing(promoter_radius)){promoter_radius = 200}
if(missing(promoter_radius2)){promoter_radius2 = 100}
if(missing(SNP)){SNP = "False"}
if(missing(adjacent)){adjacent= "False"}
if(missing(PromoterStop)){PromoterStop= "False"}
if(missing(NearestTwoDirection)){NearestTwoDirection= "True"}
if(missing(UTR3)){UTR3= "False"}
if(missing(collapsemethod)){collapsemethod= "MaxMean"}
if(missing(B)){B= 100}
if(missing(alpha)){alpha= 5}
if(missing(logCheck)){logCheck= FALSE}
if(missing(na.rm)){na.rm= FALSE}
if(missing(min_Intersect_Count)){min_Intersect_Count=5}

if(SNP %in% c("T","TRUE","True",TRUE)){SNP= "True"}
if(SNP %in% c("F","FALSE","False",FALSE)){SNP= "False"}

if(PromoterStop %in% c("T","TRUE","True",TRUE)){PromoterStop= "True"}
if(PromoterStop %in% c("F","FALSE","False",FALSE)){PromoterStop= "False"}

if(NearestTwoDirection %in% c("T","TRUE","True",TRUE)){NearestTwoDirection= "True"}
if(NearestTwoDirection %in% c("F","FALSE","False",FALSE)){NearestTwoDirection= "False"}

if(UTR3 %in% c("T","TRUE","True",TRUE)){UTR3= "True"}
if(UTR3 %in% c("F","FALSE","False",FALSE)){UTR3= "False"}
    
if(adjacent %in% c("T","TRUE","True",TRUE)){adjacent= "True"}
if(adjacent %in% c("F","FALSE","False",FALSE)){adjacent= "False"}
if(adjacent== "True"){search_radius =0}

if(length(genome>1)){genome=genome[1]}

if(!collapsemethod %in% c("MaxMean","function","ME","maxRowVariance","MinMean","absMinMean","absMaxMean","Average")){stop("please check the collapsemethod")}

data(GO_BP_list,package="seq2pathway.data")
data(GO_MF_list,package="seq2pathway.data")
data(GO_CC_list,package="seq2pathway.data")
data(Des_BP_list,package="seq2pathway.data")
data(Des_CC_list,package="seq2pathway.data")
data(Des_MF_list,package="seq2pathway.data")

#######seq2gene function
seq2gene_result<-runseq2gene(inputfile=inputfile,
            search_radius=search_radius,promoter_radius=promoter_radius,promoter_radius2=promoter_radius2,genome=genome,
            adjacent=adjacent,SNP=SNP,PromoterStop=PromoterStop,NearestTwoDirection=NearestTwoDirection,UTR3=UTR3)

seq2gene_result_fornext<-seq2gene_result[[2]]
seq2gene_result_fornext<-seq2gene_result_fornext[,c(1,13)]

genename<-unique(seq2gene_result_fornext[,2])

print("Start test..............")
############Fisher test
if(FisherTest==TRUE){
if(DataBase=="GOterm"){
FS_test<-FisherTest_GO_BP_MF_CC(gs=as.vector(genename),genome=genome,min_Intersect_Count=min_Intersect_Count)
}else{
FS_test<-FisherTest_MsigDB(gsmap=DataBase,gs=as.vector(genename),genome=genome,min_Intersect_Count=min_Intersect_Count)
}

#print("Fisher's exact test is done")
}


#############################rungene2pathway,normalization,empiricalP,summary table
if(FAIMETest==TRUE){
#Preprocess peak score with gene
    tmpinfile = tempfile()
    if(class(inputfile)=="data.frame"){
      if(ncol(inputfile)<5){stop("please check the format of the input data, some column is missing")}
      write.table(inputfile, file=tmpinfile,sep="\t",quote=FALSE,row.names = FALSE)}
    if(class(inputfile)=="GRanges"){
      test<-as.data.frame(inputfile)
      if(ncol(test)<7){stop("please check the format of the input data, some column is missing")}
      if(ncol(test)>=7){write.table(test[,c(6,1,2,3,7)], file=tmpinfile,sep="\t",quote=FALSE,row.names = FALSE)}
      }
      
peak_fornext<-read.table(file=tmpinfile,header=TRUE,sep="\t")
if(ncol(peak_fornext)<5){stop("Please check input file format, some required column is missing.")}
peak_fornext<-peak_fornext[,c(1,5)]
peak_anno_score<-merge(seq2gene_result_fornext,peak_fornext,by=names(seq2gene_result_fornext)[1],all=TRUE)
############collapse function
dat_collapsed<-Peak_Gene_Collapse(input=peak_anno_score,collapsemethod=collapsemethod)
dat_CP<-data.frame(dat_collapsed[,c(2:ncol(dat_collapsed))])
rownames(dat_CP)<-rownames(dat_collapsed)
colnames(dat_CP)<-colnames(dat_collapsed)[2:ncol(dat_collapsed)]

dat_collapsed$gene<-as.vector(toupper(rownames(dat_collapsed)))
if(DataBase=="GOterm"){
#FAIME
GO_MF_FAIME<-rungene2pathway(dat=dat_CP,gsmap=GO_MF_list,alpha=alpha,logCheck=logCheck,method="FAIME",na.rm=na.rm)
GO_BP_FAIME<-rungene2pathway(dat=dat_CP,gsmap=GO_BP_list,alpha=alpha,logCheck=logCheck,method="FAIME",na.rm=na.rm)
GO_CC_FAIME<-rungene2pathway(dat=dat_CP,gsmap=GO_CC_list,alpha=alpha,logCheck=logCheck,method="FAIME",na.rm=na.rm)
#Normalization
GO_MF_FAIME_N<-Normalize_F(input=GO_MF_FAIME)
GO_BP_FAIME_N<-Normalize_F(input=GO_BP_FAIME)
GO_CC_FAIME_N<-Normalize_F(input=GO_CC_FAIME)
#EmpiricalP
GO_MF_FAIME_Pvalue<-rungene2pathway_EmpiricalP(dat=dat_CP,gsmap=GO_MF_list,alpha=alpha,logCheck=logCheck,method="FAIME",B=B,na.rm=na.rm)
GO_BP_FAIME_Pvalue<-rungene2pathway_EmpiricalP(dat=dat_CP,gsmap=GO_BP_list,alpha=alpha,logCheck=logCheck,method="FAIME",B=B,na.rm=na.rm)
GO_CC_FAIME_Pvalue<-rungene2pathway_EmpiricalP(dat=dat_CP,gsmap=GO_CC_list,alpha=alpha,logCheck=logCheck,method="FAIME",B=B,na.rm=na.rm)

########gene2pathway table
GO_MF_N_P<-merge(GO_MF_FAIME_N,GO_MF_FAIME_Pvalue,by="row.names",all=TRUE)
rownames(GO_MF_N_P)<-GO_MF_N_P$Row.names
GO_MF_N_P<-GO_MF_N_P[,-1]
for(i in 1:nrow(GO_MF_N_P)){
intsect<-intersect(toupper(rownames(dat_CP)),toupper(unlist(GO_MF_list[names(GO_MF_list)==rownames(GO_MF_N_P)[i]])))
GO_MF_N_P$Intersect_Count[i]<-length(intsect)
GO_MF_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
GO_MF_N_P$Intersect_element[i]<-paste(as.vector(dat_collapsed[dat_collapsed$gene %in% intsect,c(1)]),collapse=" ")
rm(intsect)
GO_MF_N_P$Des[i]<-as.character(Des_MF_list[which(names(Des_MF_list)==rownames(GO_MF_N_P)[i])])
}
GO_MF_N_P<-GO_MF_N_P[GO_MF_N_P$Intersect_Count>=min_Intersect_Count,]
GO_MF_N_P<-GO_MF_N_P[,c(ncol(GO_MF_N_P),1:(ncol(GO_MF_N_P)-1))]

GO_BP_N_P<-merge(GO_BP_FAIME_N,GO_BP_FAIME_Pvalue,by="row.names",all=TRUE)
rownames(GO_BP_N_P)<-GO_BP_N_P$Row.names
GO_BP_N_P<-GO_BP_N_P[,-1]
for(i in 1:nrow(GO_BP_N_P)){
intsect<-intersect(toupper(rownames(dat_CP)),toupper(unlist(GO_BP_list[names(GO_BP_list)==rownames(GO_BP_N_P)[i]])))
GO_BP_N_P$Intersect_Count[i]<-length(intsect)
GO_BP_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
GO_BP_N_P$Intersect_element[i]<-paste(as.vector(dat_collapsed[dat_collapsed$gene %in% intsect,c(1)]),collapse=" ")
rm(intsect)
GO_BP_N_P$Des[i]<-as.character(Des_BP_list[which(names(Des_BP_list)==rownames(GO_BP_N_P)[i])])
}
GO_BP_N_P<-GO_BP_N_P[GO_BP_N_P$Intersect_Count>=min_Intersect_Count,]
GO_BP_N_P<-GO_BP_N_P[,c(ncol(GO_BP_N_P),1:(ncol(GO_BP_N_P)-1))]

GO_CC_N_P<-merge(GO_CC_FAIME_N,GO_CC_FAIME_Pvalue,by="row.names",all=TRUE)
rownames(GO_CC_N_P)<-GO_CC_N_P$Row.names
GO_CC_N_P<-GO_CC_N_P[,-1]
for(i in 1:nrow(GO_CC_N_P)){
intsect<-intersect(toupper(rownames(dat_CP)),toupper(unlist(GO_CC_list[names(GO_CC_list)==rownames(GO_CC_N_P)[i]])))
GO_CC_N_P$Intersect_Count[i]<-length(intsect)
GO_CC_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
GO_CC_N_P$Intersect_element[i]<-paste(as.vector(dat_collapsed[dat_collapsed$gene %in% intsect,c(1)]),collapse=" ")
rm(intsect)
GO_CC_N_P$Des[i]<-as.character(Des_CC_list[which(names(Des_CC_list)==rownames(GO_CC_N_P)[i])])
}
GO_CC_N_P<-GO_CC_N_P[GO_CC_N_P$Intersect_Count>=min_Intersect_Count,]
GO_CC_N_P<-GO_CC_N_P[,c(ncol(GO_CC_N_P),1:(ncol(GO_CC_N_P)-1))]

gene2pathway_result<-list()
gene2pathway_result[[1]]<-GO_BP_N_P
names(gene2pathway_result)[1]<-c("GO_BP")
gene2pathway_result[[2]]<-GO_CC_N_P
names(gene2pathway_result)[2]<-c("GO_CC")
gene2pathway_result[[3]]<-GO_MF_N_P
names(gene2pathway_result)[3]<-c("GO_MF")

}else{
dat_FAIME<-rungene2pathway(dat=dat_CP,gsmap=DataBase,alpha=alpha,logCheck=logCheck,method="FAIME",na.rm=na.rm)
N_FAIME<-Normalize_F(input=dat_FAIME)
dat_FAIME_Pvalue<-rungene2pathway_EmpiricalP(dat=dat_CP,gsmap=DataBase,alpha=alpha,logCheck=logCheck,method="FAIME",B=B,na.rm=na.rm)
N_FAIME<-as.matrix(N_FAIME)
dat_FAIME_Pvalue<-as.matrix(dat_FAIME_Pvalue)
DB_N_P<-cbind(N_FAIME,as.matrix(dat_FAIME_Pvalue[match(rownames(N_FAIME), rownames(dat_FAIME_Pvalue)),]))
colnames(DB_N_P)<-c("score2pathscore_Normalized","score2pathscore_Pvalue")
DB_N_P<-as.data.frame(DB_N_P)

for(i in 1:nrow(DB_N_P)){
if(class(DataBase)=="GSA.genesets"){
intsect<-intersect(toupper(rownames(dat_CP)),toupper(unlist(DataBase$genesets[which(DataBase$geneset.names==rownames(DB_N_P)[i])])))
}else if(class(DataBase)=="list"){
intsect<-intersect(toupper(rownames(dat_CP)),toupper(unlist(DataBase[names(DataBase)==rownames(DB_N_P)[i]])))
}
DB_N_P$Intersect_Count[i]<-length(intsect)
DB_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
DB_N_P$Intersect_element[i]<-paste(as.vector(dat_collapsed[dat_collapsed$gene %in% intsect,c(1)]),collapse=" ")
rm(intsect)
if(class(DataBase)=="GSA.genesets"){
DB_N_P$Des[i]<-as.character(DataBase$geneset.descriptions[which(DataBase$geneset.names==rownames(DB_N_P)[i])])}
}
DB_N_P<-DB_N_P[DB_N_P$Intersect_Count>=min_Intersect_Count,]
gene2pathway_result<-DB_N_P[,c(ncol(DB_N_P),1:(ncol(DB_N_P)-1))]
}
print("gene2pathway analysis is done")
}######ends:FAIMETest==TRUE

if(exists("gene2pathway_result")&exists("FS_test")){
TotalResult<-list()
TotalResult[[1]]<-seq2gene_result
names(TotalResult)[1]<-"seq2gene_result"
TotalResult[[2]]<-gene2pathway_result
names(TotalResult)[2]<-"gene2pathway_result.FAIME"
TotalResult[[3]]<-FS_test
names(TotalResult)[3]<-"gene2pathway_result.FET"
TotalResult[[4]]<-dat_CP
names(TotalResult)[4]<-"gene_collapse"
}else if(exists("gene2pathway_result")&exists("FS_test")==FALSE){
TotalResult<-list()
TotalResult[[1]]<-seq2gene_result
names(TotalResult)[1]<-"seq2gene_result"
TotalResult[[2]]<-gene2pathway_result
names(TotalResult)[2]<-"gene2pathway_result.FAIME"
TotalResult[[3]]<-dat_CP
names(TotalResult)[3]<-"gene_collapse"
}
else if(exists("gene2pathway_result")==FALSE&exists("FS_test")){
TotalResult<-list()
TotalResult[[1]]<-seq2gene_result
names(TotalResult)[1]<-"seq2gene_result"
TotalResult[[2]]<-FS_test
names(TotalResult)[2]<-"gene2pathway_result.FET"
}
return(TotalResult)
}


#function 16
gene2pathway_test<-function(dat,DataBase="GOterm",FisherTest=TRUE,EmpiricalTest=FALSE,method=c("FAIME","KS-rank","cumulative-rank"),genome=c("hg38","hg19","mm10","mm9"),alpha=5, logCheck=FALSE,na.rm=FALSE,B=100,min_Intersect_Count=5){
options(warn=-1)
if(missing(FisherTest)){FisherTest=TRUE}
if(missing(DataBase)){DataBase="GOterm"}
if(missing(genome)){genome= "hg19"}
if(length(genome>1)){genome=genome[1]}
if(missing(method)){method="FAIME"}
if(missing(alpha)){alpha=5}
if(missing(na.rm)){na.rm=FALSE}
if(missing(logCheck)){logCheck=FALSE}
if(missing(EmpiricalTest)){EmpiricalTest=FALSE}
if(missing(B)){B=100}
if(missing(min_Intersect_Count)){min_Intersect_Count=5}

data(GO_BP_list,package="seq2pathway.data")
data(GO_MF_list,package="seq2pathway.data")
data(GO_CC_list,package="seq2pathway.data")
data(Des_BP_list,package="seq2pathway.data")
data(Des_CC_list,package="seq2pathway.data")
data(Des_MF_list,package="seq2pathway.data")
############Fisher test
if(FisherTest==TRUE){
if(DataBase=="GOterm"){
FS_test<-FisherTest_GO_BP_MF_CC(gs=as.vector(rownames(dat)),genome=genome,min_Intersect_Count=min_Intersect_Count)
}else{
FS_test<-FisherTest_MsigDB(gsmap=DataBase,gs=as.vector(rownames(dat)),genome=genome,min_Intersect_Count=min_Intersect_Count)
}
#print("Fisher's exact test is done")
}

#############################rungene2pathway,normalization,empiricalP,summary table
if(DataBase=="GOterm"){
#method
GO_MF_method<-rungene2pathway(dat=dat,gsmap=GO_MF_list,alpha=alpha,logCheck=logCheck,method=method,na.rm=na.rm)
GO_BP_method<-rungene2pathway(dat=dat,gsmap=GO_BP_list,alpha=alpha,logCheck=logCheck,method=method,na.rm=na.rm)
GO_CC_method<-rungene2pathway(dat=dat,gsmap=GO_CC_list,alpha=alpha,logCheck=logCheck,method=method,na.rm=na.rm)
#Normalization
GO_MF_method_N<-Normalize_F(input=GO_MF_method)
GO_BP_method_N<-Normalize_F(input=GO_BP_method)
GO_CC_method_N<-Normalize_F(input=GO_CC_method)
#EmpiricalP
if(EmpiricalTest==TRUE){
GO_MF_method_Pvalue<-rungene2pathway_EmpiricalP(dat=dat,gsmap=GO_MF_list,alpha=alpha,logCheck=logCheck,method=method,B=B,na.rm=na.rm)
GO_BP_method_Pvalue<-rungene2pathway_EmpiricalP(dat=dat,gsmap=GO_BP_list,alpha=alpha,logCheck=logCheck,method=method,B=B,na.rm=na.rm)
GO_CC_method_Pvalue<-rungene2pathway_EmpiricalP(dat=dat,gsmap=GO_CC_list,alpha=alpha,logCheck=logCheck,method=method,B=B,na.rm=na.rm)
GO_MF_N_P<-merge(GO_MF_method_N,GO_MF_method_Pvalue,by="row.names",all=TRUE)
rownames(GO_MF_N_P)<-GO_MF_N_P$Row.names
GO_MF_N_P<-GO_MF_N_P[,-1]
GO_BP_N_P<-merge(GO_BP_method_N,GO_BP_method_Pvalue,by="row.names",all=TRUE)
rownames(GO_BP_N_P)<-GO_BP_N_P$Row.names
GO_BP_N_P<-GO_BP_N_P[,-1]
GO_CC_N_P<-merge(GO_CC_method_N,GO_CC_method_Pvalue,by="row.names",all=TRUE)
rownames(GO_CC_N_P)<-GO_CC_N_P$Row.names
GO_CC_N_P<-GO_CC_N_P[,-1]
}else{
GO_MF_N_P<-GO_MF_method_N
GO_BP_N_P<-GO_BP_method_N
GO_CC_N_P<-GO_CC_method_N
}

########gene2pathway table
for(i in 1:nrow(GO_MF_N_P)){
intsect<-intersect(toupper(rownames(dat)),toupper(unlist(GO_MF_list[names(GO_MF_list)==rownames(GO_MF_N_P)[i]])))
GO_MF_N_P$Intersect_Count[i]<-length(intsect)
GO_MF_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
rm(intsect)
GO_MF_N_P$Des[i]<-as.character(Des_MF_list[which(names(Des_MF_list)==rownames(GO_MF_N_P)[i])])
}
GO_MF_N_P<-GO_MF_N_P[GO_MF_N_P$Intersect_Count>=min_Intersect_Count,]
GO_MF_N_P<-GO_MF_N_P[,c(ncol(GO_MF_N_P),1:(ncol(GO_MF_N_P)-1))]

for(i in 1:nrow(GO_BP_N_P)){
intsect<-intersect(toupper(rownames(dat)),toupper(unlist(GO_BP_list[names(GO_BP_list)==rownames(GO_BP_N_P)[i]])))
GO_BP_N_P$Intersect_Count[i]<-length(intsect)
GO_BP_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
rm(intsect)
GO_BP_N_P$Des[i]<-as.character(Des_BP_list[which(names(Des_BP_list)==rownames(GO_BP_N_P)[i])])
}
GO_BP_N_P<-GO_BP_N_P[GO_BP_N_P$Intersect_Count>=min_Intersect_Count,]
GO_BP_N_P<-GO_BP_N_P[,c(ncol(GO_BP_N_P),1:(ncol(GO_BP_N_P)-1))]

for(i in 1:nrow(GO_CC_N_P)){
intsect<-intersect(toupper(rownames(dat)),toupper(unlist(GO_CC_list[names(GO_CC_list)==rownames(GO_CC_N_P)[i]])))
GO_CC_N_P$Intersect_Count[i]<-length(intsect)
GO_CC_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
rm(intsect)
GO_CC_N_P$Des[i]<-as.character(Des_CC_list[which(names(Des_CC_list)==rownames(GO_CC_N_P)[i])])
}
GO_CC_N_P<-GO_CC_N_P[GO_CC_N_P$Intersect_Count>=min_Intersect_Count,]
GO_CC_N_P<-GO_CC_N_P[,c(ncol(GO_CC_N_P),1:(ncol(GO_CC_N_P)-1))]

gene2pathway_result<-list()
gene2pathway_result[[1]]<-GO_BP_N_P
names(gene2pathway_result)[1]<-c("GO_BP")
gene2pathway_result[[2]]<-GO_CC_N_P
names(gene2pathway_result)[2]<-c("GO_CC")
gene2pathway_result[[3]]<-GO_MF_N_P
names(gene2pathway_result)[3]<-c("GO_MF")

}else{
dat_method<-rungene2pathway(dat=dat,gsmap=DataBase,alpha=alpha,logCheck=logCheck,method=method,na.rm=na.rm)
N_method<-Normalize_F(input=dat_method)
if(EmpiricalTest==TRUE){
dat_method_Pvalue<-rungene2pathway_EmpiricalP(dat=dat,gsmap=DataBase,alpha=alpha,logCheck=logCheck,method=method,B=B,na.rm=na.rm)
N_method<-as.matrix(N_method)
dat_method_Pvalue<-as.matrix(dat_method_Pvalue)
if(ncol(dat_method_Pvalue)>1){
DB_N_P<-cbind(N_method,dat_method_Pvalue[match(rownames(N_method), rownames(dat_method_Pvalue)),]) }
if(ncol(dat_method_Pvalue)==1){
DB_N_P<-cbind(N_method,as.matrix(dat_method_Pvalue[match(rownames(N_method), rownames(dat_method_Pvalue)),])) 
colnames(DB_N_P)<-c("score2pathscore_Normalized","score2pathscore_Pvalue")
}
DB_N_P<-as.data.frame(DB_N_P)
}else{
DB_N_P=N_method}

for(i in 1:nrow(DB_N_P)){
if(class(DataBase)=="GSA.genesets"){
intsect<-intersect(toupper(rownames(dat)),toupper(unlist(DataBase$genesets[which(DataBase$geneset.names==rownames(DB_N_P)[i])])))
}else if(class(DataBase)=="list"){
intsect<-intersect(toupper(rownames(dat)),toupper(unlist(DataBase[names(DataBase)==rownames(DB_N_P)[i]])))
}
DB_N_P$Intersect_Count[i]<-length(intsect)
DB_N_P$Intersect_gene[i]<-paste(intsect,collapse=" ")
rm(intsect)
if(class(DataBase)=="GSA.genesets"){
DB_N_P$Des[i]<-as.character(DataBase$geneset.descriptions[which(DataBase$geneset.names==rownames(DB_N_P)[i])])}
}
DB_N_P<-DB_N_P[DB_N_P$Intersect_Count>=min_Intersect_Count,]
gene2pathway_result<-DB_N_P[,c(ncol(DB_N_P),1:(ncol(DB_N_P)-1))]
}
print("gene2pathway analysis is done")


if(exists("gene2pathway_result")&exists("FS_test")){
TResult<-list()
TResult[[1]]<-gene2pathway_result
names(TResult)[1]<-"gene2pathway_result.2"
TResult[[2]]<-FS_test
names(TResult)[2]<-"gene2pathway_result.FET"
}else if(exists("gene2pathway_result")&exists("FS_test")==FALSE){
TResult<-gene2pathway_result
}
else if(exists("gene2pathway_result")==FALSE&exists("FS_test")){
TResult<-FS_test
}
return(TResult)

}


