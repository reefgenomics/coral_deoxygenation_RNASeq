setwd()
library(DESeq2)
library(ggplot2)
library(ellipse)

######### pca plots ########
map.pre<-read.delim("./Input_files/Metadata.txt", row.names = 1, header = TRUE)
map<-map.pre[order(rownames(map.pre)), ]
cnt<-read.table("./Input_files/Finaltable_species", header=TRUE, row.names = 1, sep= "\t")
cnts<-round(cnt, 0)
cnts.sub<-subset(cnts, rowSums(cnts) > 50) 
cnts.sor<-cnts.sub[ , order(names(cnts.sub))]
colnames(cnts.sor)=gsub("\\.", "-", colnames(cnts.sor))

map$Group <- paste(map$Condition, map$Species, sep="_")
map$Group <- paste(map$Condition, map$Species, map$Time, sep="_")
map.a<-subset(map, map$Time =="T3")
cnt.a<-cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.a))])] #make names to match up in same order
dds.a<-DESeqDataSetFromMatrix(countData = as.matrix(cnt.a), colData = map.a, design = ~Group)#Species + Condition)
dds.a<-estimateSizeFactors(dds.a)
dds.a<-estimateDispersions(dds.a)
vsd.a<-vst(dds.a)

COL2<-c("#FF99CC", "#CC0033", "#99CCFF", "#333399")
pca.vsd<-plotPCA(vsd.a, intgroup=c("Species", "Condition", "Time"), ntop = 1000, returnData = TRUE) #intgroup=interesting groups is a character vector of names in colData(x) to use for grouping
#ntop=to change how many of the most variable genes are included
#returnData= allows pca plotting with Deseq
percentVar <- round(100 * attr(pca.vsd, "percentVar"))

centroids <- aggregate(cbind(PC1,PC2)~group,pca.vsd,mean)
conf.rgn  <- do.call(rbind,lapply(unique(pca.vsd$group),function(t)
  data.frame(group=as.character(t),
             ellipse(cov(pca.vsd[pca.vsd$group==t,1:2]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.95),
             stringsAsFactors=FALSE)))
#level= to change % confidence level

ggplot(data=pca.vsd,(aes(x=PC1,y=PC2,colour = group)))+
  geom_point(size=5)+
  geom_path(data=conf.rgn, size=1.1)+
  theme_classic()+
  scale_colour_manual(values=COL2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(axis.text = element_text(size=18))+
  theme(axis.title = element_text(size=22))+
  theme(text = element_text(size=22))+
  theme(aspect.ratio = 1)

