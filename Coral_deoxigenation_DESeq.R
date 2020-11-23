
library(DESeq2)
library(ggplot2)

map.pre<-read.delim("./Input_files/Metadata.txt", row.names = 1, header = TRUE)
map<-map.pre[order(rownames(map.pre)), ]

### for condition comparison for each species ###
map.A<-subset(map, map$Species =="A.selago") #change species insert name accordingly
map.a<-map.A[order(rownames(map.A)), ] 
map1<-subset(map.a, map.a$Time == "T1") 
map3<-subset(map.a, map.a$Time == "T2") 
map5<-subset(map.a, map.a$Time == "T3") 

cnt_S<-read.table("./Input_files/Finaltable_species", header=TRUE, row.names = 1, sep="\t") 
#remove genes uniquely mapped to reference by one species
rm= read.table("./Input_files/DiffMapGenes_1433_toRemove.txt", header=TRUE, row.names = 1, sep= "\t")
cnt_2=subset(cnt_S, !(row.names(cnt_S) %in% row.names(rm)))
cnts=round(cnt_2, 0)
cnts.sub<-subset(cnts, rowSums(cnts) > 50) # 10 excludes 200k, 20 excludes almost 500k
cnts.sor<-cnts.sub[ , order(names(cnts.sub))]
colnames(cnts.sor)=gsub("\\.", "-", colnames(cnts.sor))
cnt1<-cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map1))])]
cnt3<-cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map3))])]
cnt5<-cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map5))])]

#comparison 1
dds1<-DESeqDataSetFromMatrix(countData = as.matrix(cnt1), colData = map1, design = ~Condition)
dds1<- DESeq(dds1)
res1<-results(dds1, contrast=c("Condition","treatment","control"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res1, alpha=0.05), "comparison 1: T1, treatment vs control")
results1<-as.data.frame(subset(res1, res1$padj<0.05))
write.table(results1, "output1.txt", sep = "\t", quote = FALSE, row.names= TRUE)


#comparison 3 
dds3<-DESeqDataSetFromMatrix(countData = as.matrix(cnt3), colData = map3, design = ~Condition)
dds3<-DESeq(dds3)
res3<-results(dds3, contrast=c("Condition","treatment","control"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res3, alpha=0.05), "comparison 3: T2, treatment vs control")
results3<-as.data.frame(subset(res3, res3$padj<0.05))
write.table(results3, "output2.txt", sep = "\t", quote = FALSE, row.names= TRUE)

#comparison 5
dds5<-DESeqDataSetFromMatrix(countData = as.matrix(cnt5), colData = map5, design = ~Condition)
dds5<- DESeq(dds5)
res5<-results(dds5, contrast=c("Condition","treatment","control"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res5, alpha=0.05), "comparison 5: T3, treatment vs control")
results5<-as.data.frame(subset(res5, res5$padj<0.05))
write.table(results5, "output3.txt", sep = "\t", quote = FALSE, row.names= TRUE)

