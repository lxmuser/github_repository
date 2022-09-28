#final Rscript for CA1
#team member: Li Xinmiao(A0251308M),Zhang Jingyi(A0251345L)
#data we use: GSE71348 (although it is listed in human group in CA1 pdf, actually it belongs to mouse species )

library(GEOquery)
library(limma)
library(oligo)
library(pheatmap)


#####obtain data

gse71348 <- getGEO('gse71348',GSEMatrix = TRUE,getGPL = FALSE,destdir = "./") #45101 features,15 samples

#get to know this data
class(gse71348) #list
View(gse71348) #1 list
gse71348[[1]] #extract gse71348. 45101 probe, 15 sample, GPL1261-56135
class(gse71348[[1]]) #expressionset
expressdata <- exprs(gse71348[[1]]) #get sample expression matrix
pd <- pData(gse71348[[1]]) #get sample clinical information，eg.treatment, cell line, etc.

#make it more accessible
names(pd)
names(pd)[39] <- 'treatment' #2i rspo3 chir
names(pd)[38] <- 'timepoint' #d0 d4 d6
pd$group <- as.factor(paste(pd$treatment,pd$timepoint))
levels(pd$group) #look the order
levels(pd$group) <- c('control','chir_4d','chir_6d','rspo3_4d','rspo3_6d') #rename levels in known order
grouplist <- c("control","control","control",
               "rspo3_4d","rspo3_4d","rspo3_4d",
               "chir_4d","chir_4d","chir_4d",
               "rspo3_6d","rspo3_6d","rspo3_6d",
               "chir_6d","chir_6d","chir_6d") #because there are many groups and few samples, it is not as convenient to use ifelse as it is to use this directly.
save(expressdata,grouplist,file = 'expreset.Rdata') #it can be constructed to facilitate subsequent direct extraction of data for processing if this data is well normalised and background noise does not need to be removed.


#####check the availability of data

#check the distribution of exprs(gse71348)
options(stringsAsFactors = F) #false for string as factor
load(file='expreset.Rdata') 
boxplot(expressdata,las=2,col="red") #data is uneven，need do background correction and normalization


#####pre-processing of the raw data: background correction and normalization

#get raw data:Manually download the original GSE71348_RAW.tar file from the GEO database and extract it. Also, download the probe information file GPL1261.txt for the chip used in the GSE71348. Place the above files in the same folder. 
list.files("d:/R/w6/bookcamp/gse71348_affy_expression/") #.txt is the annotation file. .cel is the raw data of sample processed by chip
celfile <- list.celfiles("d:/R/w6/bookcamp/gse71348_affy_expression/",full.name=T) #put cel file together
raw.cel <- oligo::read.celfiles(celfile) #read in cel file

#get to know this raw data
image(raw.cel) #since raw data is the signal intensity, so we can get to know this raw data by imaging it directly
rawdata <- raw.cel
exprs(rawdata)[1:10,1:10] #see the intensity of the photoelectric before precessing(raw data),every column is a sample
max(exprs(rawdata)) #see the biggest intensity is 65534. Since 256*256=65536 and two number can't be shown,so the biggest intensity is 65534 which also indicates us the raw data need to be pre-processed
pData(rawdata) #to facilitate DESeq analysis, construct pData
pData(rawdata)$filename <- sampleNames(rawdata) #save in filename (samplename)
pData(rawdata)$group <- c("control","control","control","rspo3_4d","rspo3_4d","rspo3_4d","chir_4d","chir_4d","chir_4d","rspo3_6d","rspo3_6d","rspo3_6d","chir_6d","chir_6d","chir_6d") #add group name

#visualize raw data before processing
oligo::hist(rawdata,lwd=2,xlab='log intensity',main="CEL file densities before background correcting and normalization")
boxplot(rawdata,las=2,col="red")

#visualize raw data after processing
normdata <- oligo::rma(rawdata,background=T,normalize=T) #background correction performed within the sample and normalization performed between samples.annotation:pd.mouse430.2 
oligo::hist(normdata,lwd=2,xlab='log intensity',main="CELl file densities after background correcting and normalization")
boxplot(normdata,las=2,col="red") # the medium is corrected to be at approximately  the same level

#save new data
expressdata1 <- exprs(normdata) #get new sample expression matrix
colnames(expressdata1) <- colnames(expressdata)
pd1 <- pData(normdata) #get new clinical information
save(expressdata1,grouplist,file = 'expreset1.Rdata')


#####principle component analysis(PCA) and hierarchical aggregation analysis(HAA)

#PCA
load(file = 'expreset1.Rdata')
expressdata1[1:4,1:4] #check data
expressdata1t <- t(expressdata1) #drawing PCA requeires rows to be sample names and columns to be probe names, but the expressdata1 now is opposit, so transposition is needed.
expressdata1t[1:4,1:4] #has been transposed
expressdata1t <- as.data.frame(expressdata1t) #convert matrix to dataframe
expressdata1t <- cbind(expressdata1t,grouplist) #add grouplist to the last column--is not useful in PCA but helpful to recognize sample
install.packages("FactoMineR")
install.packages("factoextra") #2 packages needed to draw PCA
library("FactoMineR")
library("factoextra")
dat.pca <- PCA(expressdata1t[,-ncol(expressdata1t)],graph = FALSE) #doing PCA without the last column
head(dat.pca$ind$coord) #numbers of probes is reduced to 5 dimensions
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = expressdata1t$grouplist,
             addEllipses = TRUE,
             legend.title="groups")
#but we think is not useful for our data.The reason is that PCA is for comparing and classifying probes and making probes dimensionality reduction. After doing PCA, we only have probes in 5 dimensionalities with 15 samples. Thus, it will be more suitable for cases with many parallel samples.
#expressdata1 has not been manipulated

#HAA
load(file = 'expreset1.Rdata')
expressdata1[1:4,1:4] #每次都要检查数据 
expressdata1t <- t(expressdata1) #drawing PCA requeires rows to be probe names and columns to be sample names, so expressdata1 is what we need.
d <- dist(expressdata1t) #calculating the distance of samples with similar gene expression
fit.complete <- hclust(d,method = "complete") #doing cluster
plot(fit.complete,hang = -1,cex=0.8) #hanging labels over the axes
#In HAA, two samples with similar expressions will be clustered, then compared with others and clustered again until all the samples are clustered. As our data is a small sample size and the number of repetitions is also small, it can be grouped manually instead of using HAA.


#####annotation for probe ID
BiocManager::install("mouse430a2.db") #Download the annotation package mouse430a2.db which corresponds to the probe ID in platform GPL1261
library(mouse430a2.db)
ids <- toTable(mouse430a2SYMBOL) #extraction of genes corresponding to the probe ID of the GPL1261 chip platform. 21402
expressdata1 <-as.data.frame(expressdata1)
expressdata1$probe_id <- rownames(expressdata1) #add probe ID as a new column
expressdata1 <- merge(expressdata1,ids,by='probe_id') #Matching and combining the probe ID in the annotation file with the probe ID in the expression matrix
#after matching, the number of features changed from 45101 to 21402 which means a significant proportion of probes did not have corresponding gene symbols and were automatically ignored, leaving only those probes with corresponding gene names.

expressdata1 <- avereps(expressdata1,ID=expressdata1$symbol) #Averaging the gene expression values for the same gene calculated by multiple probes and combining them into one column.21402-->13091
expressdata1 <- as.data.frame(expressdata1) 
row.names(expressdata1) <- expressdata1$symbol #changing line names to gene symbols instead of 1,2,3,4
expressdata1 <- expressdata1[,-c(1,17)] #the first column is the probe ID and the last column is the symbol, not the gene symbol.So remove the first and last column and keep only the gene expression.

expressdata_1 <- data.frame(expressdata1,stringsAsFactors = F) #what was a character is still a charater
expressdata_1 <- data.matrix(expressdata_1) #change frameto matrix
save(expressdata_1,grouplist,file = 'expreset2.Rdata')


#####differential expression analysis and visualization

#exam the data
load(file = 'expreset2.Rdata')
table(grouplist) 
boxplot(expressdata_1[1,]~grouplist) #boxplot can be drawn for each gene expression by grouping
t.test(expressdata_1[1,]~grouplist) #more than two groups of data are not applicable to t-test and simply design

#differential expression analysis: create contrasts matrix
design <- model.matrix(~ 0 + pd$group) 
colnames(design) <- levels(pd$group)
contrasts_matrix <- makeContrasts(rspo3_in_es_4d = rspo3_4d - control,
                                  rspo3_in_es_6d = rspo3_6d - control,
                                  chir_in_es_4d = chir_4d - control,
                                  chir_in_es_6d = chir_6d - control,
                                  interaction = (rspo3_6d - chir_6d) - (rspo3_4d - chir_4d),
                                  levels = design) 
fit <- lmFit(expressdata_1,design) #lmFitfits a linear model model given a eries of arrays for each gene
fit2 <- contrasts.fit(fit,contrasts = contrasts_matrix) 
fit3 <- eBayes(fit2) #eBayes gives a microarray linear model fit, calculated by Bayesian adjustment of the standard error to a common value
summary(decideTests(fit3,lfc=1))
topTable(fit3)

#see the number of differential expression gene without matching probe ID which means have 45101 features
expressdata1 <- exprs(normdata) #获取样本表达矩阵。行：probeid，列：sampleid
colnames(expressdata1) <- colnames(expressdata)
expressdata1 <- data.frame(expressdata1,stringsAsFactors = F)
expressdata1 <- data.matrix(expressdata1)
fit0 <- lmFit(expressdata1,design)
fit9 <- contrasts.fit(fit0,contrasts = contrasts_matrix)
fit8 <- eBayes(fit9)
summary(decideTests(fit8,lfc=1))

#heatmap
library(RColorBrewer)
heatmap(expressdata_1,
        col=rev(brewer.pal(10,"RdBu")),
        distfun = function(x)as.dist(1-cor(t(x))),
        scale = "row") #heatmap of gene expression
heatmap(fit3,
        col=rev(brewer.pal(10,"RdBu")),
        distfun = function(x)as.dist(1-cor(t(x))),
        scale = "row") #heatmap of gene expression under contrasts

#volcano
volcanoplot(fit3,xlab='rspo3_in_es_4d')

library(ggplot2)
tempoutput <- topTable(fit3,n=Inf) 
tempoutput <- na.omit(tempoutput) #remove NA values, is data.frame
tempoutput1 <- tempoutput
tempoutput1$p <- -log10(tempoutput1$P.Value)
tempoutput1$groupr4 <- ifelse(tempoutput1$P.Value > 0.01,'stable',
                              ifelse(tempoutput1$rspo3_in_es_4d > 1.5,'up',
                                     ifelse(tempoutput1$rspo3_in_es_4d < 1.5,'down','stable')))
#groupr6,groupc4,groupc6,groupin can also be set up by changing the code of tempoutput$repo3_in_es_4d part.

ggplot(data=tempoutput1,aes(x=rspo3_in_es_4d,y=p,color=groupr4)) + geom_point()
#but it seems like heatmaps more visualizable