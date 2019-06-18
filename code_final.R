library(dplyr)
library(fossil)
library(stringr)
library(superheat)
library(viridis)
library(RColorBrewer)
library(Rtsne)
library(plyr)
library(wesanderson)



#load the data table with tpm values per sample
genexp2 <- read.table("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/genexp2.tsv")

###calculate basic statistics 
genexp2$meanexp <- rowMeans(genexp2[,2:21])
genexp2$maxexp <- apply(genexp2[,2:21],1, max)
genexp2$minexp <- apply(genexp2[,2:21],1, min)
genexp2$sd <- apply(genexp2[,2:21],1, sd)
genexp2$cv <- genexp2$sd/genexp2$meanexp

#add the pfam annotation for exploring. The Pfam annotation can be downloaded from xxx
pfam <- read.csv("~/Documents/Genomes_10feb16_Shared/Sarc4_5des17/Sarc4_Pfamscan.arqdom", header = F, sep ="\t")
colnames(pfam) <- c("gene", "pfam")
genexp2 <- merge(genexp2, pfam, by.x="gene", by.y="gene", all = T)


#filtering out non-expressed genes and saving the tables to be clustered by clust.
genexp2_filt <- subset (genexp2, meanexp > 0.5)
genexp2_avg <- data.frame(gene = genexp2_filt$gene, T12 = (genexp2_filt$BT12 + genexp2_filt$CT12)/2,  T18 = (genexp2_filt$BT18 + genexp2_filt$CT18)/2,  T24 = (genexp2_filt$BT24 + genexp2_filt$CT24)/2,  T30 = (genexp2_filt$BT30 + genexp2_filt$CT30)/2,  T36 = (genexp2_filt$BT36 + genexp2_filt$CT36)/2,  T42 = (genexp2_filt$BT42 + genexp2_filt$CT42)/2,  T48 = (genexp2_filt$BT48 + genexp2_filt$CT54)/2,  T54 = (genexp2_filt$BT54 + genexp2_filt$CT60)/2,  T60 = (genexp2_filt$BT60 + genexp2_filt$CT66)/2,  T66 = (genexp2_filt$BT66 + genexp2_filt$CT72)/2)
dir.create("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final")
dir.create("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/all")
dir.create("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/rep1")
dir.create("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/rep2")
dir.create("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/avg")
write.table(genexp2_filt[,1:21], file = "~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/all/forclust_all.tsv",sep = " ", row.names = F)
write.table(genexp2_filt[,1:11], file = "~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/rep1/forclust_r1.tsv",sep = " ", row.names = F)
write.table(genexp2_filt[,c(1,12:21)], file = "~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/rep2/forclust_r2.tsv",sep = " ", row.names = F)
write.table(genexp2_avg[,1:11], file = "~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/avg/forclust_avg.tsv",sep = " ", row.names = F)


# comparison of clustering results ----------------------------------------

###AFTER running clust and saving results in the right directories, load the cluster membership files
clust_all <- read.csv("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/clust_all/Clusters_Objects.tsv", header = T, sep = "\t" )
clust_all <- clust_all[-1,]
clust_r1 <- read.csv("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/clust_r1/Clusters_Objects.tsv", header = T, sep = "\t" )
clust_r1 <- clust_r1[-1,]
clust_r2 <- read.csv("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/clust_r2/Clusters_Objects.tsv", header = T, sep = "\t" )
clust_r2 <- clust_r2[-1,]
clust_avg <- read.csv("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/clust_avg/Clusters_Objects.tsv", header = T, sep = "\t" )
clust_avg <- clust_avg[-1,]

##rewriting the datasets and merging
clus_all <- data.frame(Genes = character(), clu_all = double())
for (i in c(1:ncol(clust_all))) {
  
  temp <- as.data.frame(clust_all[,i])
  temp <- as.data.frame(temp[!apply(temp == "", 1, all),])
  colnames(temp) <- "Genes"
  #write.table(temp, file = paste0("cluster",as.character(i)), row.names = F, col.names = F, quote = F) #uncomment for saving the gene list for GO enrichment results
  
  temp$clu_all <- rep(i,each = nrow(temp))
  clus_all <- rbind(clus_all,temp)
  
  #assign(paste0("clust", as.character(i)), temp)
}

clus_avg <- data.frame(Genes = character(), clu_avg = double())
for (i in c(1:ncol(clust_all))) {
  
  temp <- as.data.frame(clust_avg[,i])
  temp <- as.data.frame(temp[!apply(temp == "", 1, all),])
  colnames(temp) <- "Genes"
  #write.table(temp, file = paste0("cluster",as.character(i)), row.names = F, col.names = F, quote = F)
  
  temp$clu_avg <- rep(i,each = nrow(temp))
  clus_avg <- rbind(clus_avg,temp)
  
  #assign(paste0("clust", as.character(i)), temp)
}

clus_r1 <- data.frame(Genes = character(), clu_r1 = double())
for (i in c(1:ncol(clust_r1))) {
  
  temp <- as.data.frame(clust_r1[,i])
  temp <- as.data.frame(temp[!apply(temp == "", 1, all),])
  colnames(temp) <- "Genes"
  #write.table(temp, file = paste0("cluster",as.character(i)), row.names = F, col.names = F, quote = F)
  
  temp$clu_r1 <- rep(i,each = nrow(temp))
  clus_r1 <- rbind(clus_r1,temp)
  
  #assign(paste0("clust", as.character(i)), temp)
}

clus_r2 <- data.frame(Genes = character(), clu_r2 = double())
for (i in c(1:ncol(clust_r2))) {
  
  temp <- as.data.frame(clust_r2[,i])
  temp <- as.data.frame(temp[!apply(temp == "", 1, all),])
  colnames(temp) <- "Genes"
  #write.table(temp, file = paste0("cluster",as.character(i)), row.names = F, col.names = F, quote = F)
  
  temp$clu_r2 <- rep(i,each = nrow(temp))
  clus_r2 <- rbind(clus_r2,temp)
  
  #assign(paste0("clust", as.character(i)), temp)
}

##analysis of agreement of clustering methods
clus_membership <- merge(clus_r2, (merge(clus_r1,(merge (clus_all, clus_avg, by="Genes", all = T)), by="Genes", all = T)),by="Genes", all = T)
clus_membership[is.na(clus_membership)] = 0

plot(jitter(clus_membership$clu_avg), jitter(clus_membership$clu_all), pch = 16, cex = 0.5, xlab = "cluster membership (averaged replicates)", ylab = "", xaxt = "n", yaxt = "n")
axis (side = 1, at = 0:9, labels = c("nc", 1:9))
axis (side = 2, at = 0:9, labels = c("nc", 1:9))

plot(jitter(clus_membership$clu_r2), jitter(clus_membership$clu_all), pch = 16, cex = 0.5, xlab = "cluster membership (replicate 2 only)", ylab = "", xaxt = "n", yaxt = "n")
axis (side = 1, at = 0:9, labels = c("nc", 1:9))
axis (side = 2, at = 0:9, labels = c("nc", 1:9))

plot(jitter(clus_membership$clu_r1), jitter(clus_membership$clu_all), pch = 16, cex = 0.5, xlab = "cluster membership (replicate 1 only)", ylab = "cluster membership (all data)", xaxt = "n", yaxt = "n")
axis (side = 1, at = 0:10, labels = c("nc", 1:10))
axis (side = 2, at = 0:9, labels = c("nc", 1:9))

rand_allvsr1 <- adj.rand.index(clus_membership$clu_all, clus_membership$clu_r1)
rand_allvsr2 <- adj.rand.index(clus_membership$clu_all, clus_membership$clu_r2)
rand_allvsavg <- adj.rand.index(clus_membership$clu_all, clus_membership$clu_avg)


# plotting and analysis of clustering results -----------------------------

#load the NORMALIZED RNA abundance table produced by clust
genexp2_norm <- read.table("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/clustdata_final/clust_all/Processed_Data/forclust_all.tsv_processed.tsv", header = T)

genexp2_filt <- merge(genexp2_filt, clus_all, by.x="gene", by.y="Genes", all = T)
genexp2_statsonly <- genexp2_filt[,c(1,22:28)]
genexp2_norm <- merge (genexp2_norm, genexp2_statsonly, by.x="Genes", by.y="gene", all = T)
genexp2_clustonly <- subset(genexp2_norm,!is.na(genexp2_norm$clu_all))

###clustering of timepoints of the filtered dataset
d <- cor(genexp2_filt[,2:21], method="spearman")
hc <- hclust(dist(1-d))
plot(hc, show.node.label=TRUE, no.margin=TRUE,xlab="",main="",hang=-0.1, ylab = "clustering of samples")

#heatmap with cluster membership of coding genes only
genexp2_clustonlyORFonly <- subset(genexp2_clustonly,str_detect(Genes, "Sarc4_g"))
pal = (cividis(100))
superheat(as.matrix(genexp2_clustonlyORFonly[,2:21]),heat.pal=pal,membership.rows = genexp2_clustonlyORFonly$clu_all)

###separate heatmaps of rep1 and rep2
superheat(as.matrix(genexp2_clustonlyORFonly[,2:11]),heat.pal=pal,membership.rows = genexp2_clustonlyORFonly$clu_all, heat.lim=c(-3.58,3.42), heat.pal.values = c(0,0.4,0.45,0.5,0.55,0.6,0.65,1))
superheat(as.matrix(genexp2_clustonlyORFonly[,12:21]),heat.pal=pal,membership.rows = genexp2_clustonlyORFonly$clu_all, heat.lim=c(-3.58,3.42), heat.pal.values = c(0,0.4,0.45,0.5,0.55,0.6,0.65,1))

###plotting averages of within-cluster mRNA abundance (both replicates averaged)
pal2c <- c(brewer.pal(4,"Blues")[c(2:4)],brewer.pal(6,"Oranges"))
clumeans <- matrix(, nrow = length(unique(genexp2_clustonlyORFonly$clu_all)), ncol = 20)
for (i in sort(unique(genexp2_clustonlyORFonly$clu_all))) {
  bbbb <- subset(genexp2_clustonlyORFonly,clu_all==i)
  clumeans[i,] <- colMeans(bbbb[,2:21])
}
clumeans1 <- clumeans[,1:10]
clumeans2 <- clumeans[,11:20]
clumeans <- 0.5*(clumeans1 + clumeans2)

plot(1:10,clumeans[1,1:10],xlab="time (h)", ylab="norm. gene expression",type="o",pch=20, lwd = 1.5,ylim=c(-3,4),col=pal2c[1], xaxt = "n")
lines(1:10,clumeans[2,1:10],col=pal2c[2], type = "o", pch=20, lwd = 1.5)
lines(1:10,clumeans[3,1:10],col=pal2c[3], type = "o", pch=20, lwd = 1.5)
lines(1:10,clumeans[4,1:10],col=pal2c[4], type = "o", pch=20, lwd = 1.5)
lines(1:10,clumeans[5,1:10],col=pal2c[5], type = "o", pch=20, lwd = 1.5)
lines(1:10,clumeans[6,1:10],col=pal2c[6], type = "o", pch=20, lwd = 1.5)
lines(1:10,clumeans[7,1:10],col=pal2c[7], type = "o", pch=20, lwd = 1.5)
lines(1:10,clumeans[8,1:10],col=pal2c[8], type = "o", pch=20, lwd = 1.5)
lines(1:10,clumeans[9,1:10],col=pal2c[9], type = "o", pch=20, lwd = 1.5)
axis(1,at = 1:10, labels <- 6*c(2:11))
legend(1,4,col = pal2c, legend = paste("cluster", c(1:9)), lty = 1, cex = 0.6, bty = "n")

##tsne plot
tsne2_cluonly <- Rtsne(genexp2_clustonlyORFonly[2:21], dims=2, perplexity = 100, max.iter = 1000, check_duplicates = FALSE)

lab2c <- genexp2_clustonlyORFonly$clu_all
labe2c <- as.factor(lab2c)
plot(tsne2_cluonly$Y, col = pal2c[labe2c], pch=20, main="t-SNE plot", xlab = "", ylab= "")


# phylostratigraphic analysis of gene expression --------------------------

##geneages file generated using the orthogroup file and Count software
geneagetable <- read.csv("~/Documents/Sphaero_orthofinderproject/Results_Nov20/geneages.csv")

genexp2_filt <- merge(genexp2_filt,geneagetable,by.x="gene",by.y="gene",all = T)
genexp2_norm <- merge(genexp2_norm,geneagetable,by.x="Genes",by.y="gene",all = T)

genexp2_filt <- subset(genexp2_filt,meanexp > 0.5)
genexp2_norm <- subset(genexp2_norm,meanexp > 0.5)
genexp2_normORFonly <- subset(genexp2_norm,str_detect(Genes,"Sarc4_g"))
genexp2_normORFonly$age[is.na(genexp2_normORFonly$age)] = 1
genexp2_normORFonly$clu_all[is.na(genexp2_normORFonly$clu_all)] = 0

###pooling together gene age strata into only 4

genexp2_normORFonly$age2 <- mapvalues(genexp2_normORFonly$age, c(1:10),c(1,2,2,2,2,3,3,3,4,4))
ages2 <- sort(unique(genexp2_normORFonly$age2))
clusters <- sort(unique(genexp2_normORFonly$clu_all))
clusters <- clusters[-1]


##calculating enrichment of individual phylostrata in gene expression clusters
pvalmat2 <- matrix(rep(1,times=length(clusters)*length(ages2)), nrow = length(clusters), ncol = length(ages2))
enrichmat2 <- matrix(rep(1,times=length(clusters)*length(ages2)), nrow = length(clusters), ncol = length(ages2)) 

for (i in clusters) {
  for (j in ages2) {
    
    aa <- matrix(c ( sum(genexp2_normORFonly$age2==j & genexp2_normORFonly$clu_all==i), sum(genexp2_normORFonly$clu_all==i & genexp2_normORFonly$age2 !=j), sum(genexp2_normORFonly$age2==j & genexp2_normORFonly$clu_all != i), sum(genexp2_normORFonly$age2 != j & genexp2_normORFonly$clu_all != i) ), nrow = 2)
    pvalmat2[i,j] <- fisher.test(aa)$p.value
    enrichmat2[i,j] <- (aa[1,1]/(aa[2,1] + aa[1,1]) - aa[1,2]/(aa[2,2]+aa[1,2]))*100
    
  }
}

##plots with age enrichments in each cluster
pal3c <- wes_palette("Zissou1",10,type="continuous")

layout(matrix((1:10),2,5,byrow = TRUE))
for (l in 1:nrow(enrichmat2)){
  par(mar=c(2.5,2.5,2.5,2.5))
  barplot(enrichmat2[l,], col = pal3c[c(1,5,8,10)], main = paste("cluster",l), ylim = c(-30,40), ylab = "rel. enrichment (%)")
}

##pie chart with fractions of genes per stratum
ages <- sort(unique(genexp2_normORFonly$age))
numgenes <- rep(1,length(ages))
for (j in ages){
  numgenes[j] <- sum(genexp2_normORFonly$age == j)
}

numgenes2 <- rep(1,length(ages2))
for (j in ages2){
  numgenes2[j] <- sum(genexp2_normORFonly$age2 == j)
}

dev.off()
pie(numgenes, labels= c("Sarc","SarcCreo", "SarcCreoIhof", "Ichthoponida", "Ichthyosporea", "Teretosporea","Holozoa", "Opisthokonta","Unikonta","Eukaryota"),col = pal3c)
pie(numgenes2, labels= c("Sarc", "Ichthyosporea", "Opisthokonta","Eukaryota"),col = pal3c[c(1,5,8,10)])


###distributions of mean expression and cv of expression per phylostratum
boxplot(genexp2_normORFonly$cv ~ genexp2_normORFonly$age,outline = F,
        col = pal3c, xlab = "", ylab = "CV", ylim=c(0,2.2),xaxt="n" )
axis(side = 1, at = 1:10, labels= c("Sarc","SarcCreo", "SarcCreoIhof", "Ichthoponida", "Ichthyosporea", "Teretosporea","Holozoa", "Opisthokonta","Unikonta","Eukaryota"),cex.axis = 0.7,las = 2)

boxplot(genexp2_normORFonly$meanexp ~ genexp2_normORFonly$age,outline = F,
        col = pal3c, xlab = "", ylab = "mean expression levels (tpm)",xaxt="n" )
axis(side = 1, at = 1:10, labels= c("Sarc","SarcCreo", "SarcCreoIhof", "Ichthoponida", "Ichthyosporea", "Teretosporea","Holozoa", "Opisthokonta","Unikonta","Eukaryota"),cex.axis = 0.7)

###transcriptional age index
genexp2_filtORFonly <- subset(genexp2_filt,str_detect(gene,"Sarc4_g"))
genexp2_filtORFonly <- subset(genexp2_filtORFonly,!is.na(genexp2_filtORFonly$BT12))
genexp2_filtORFonly$age[is.na(genexp2_filtORFonly$age)] = 1

genexp2_phyloindex <- (genexp2_filtORFonly[,c(2:21)]*(11-genexp2_filtORFonly$age))
hourglass <- colSums(genexp2_phyloindex)/colSums(genexp2_filtORFonly[,2:21])

plot(hourglass[1:10],xlab = "time (h)",ylab = "Transcriptome age index", type ="o", xaxt = "n")
axis(1,at = 1:10, labels <- 6*c(2:11))
lines(hourglass[11:20], lty = 4)
points(hourglass[11:20], pch = 16)
legend(1,3.2, lty = c(1,4), pch =c(1,16), legend = c("replicate 1","replicate 2"),bty = "n")


# lincRNA expression analysis ---------------------------------------------

conservedLncs2 <- read.csv2("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/Conserved_LncRNA_list_Abeo_Pirum_Cfra.txt", header = F, sep = "\t")
conservedLncs <- read.csv("~/Dropbox/Sarc_timecourse_mRNAseq_newindex/Conserved_LncRNA_list.txt", header = F, sep = "\t" )
genexp2_norm_lnc_cons <- merge(genexp2_norm, conservedLncs, by.x="Genes", by.y="V1", all = F)
genexp2_norm_lnc_cons <- subset(genexp2_norm_lnc_cons)
genexp2_norm_lnc_cons2 <- merge(genexp2_norm, conservedLncs2, by.x="Genes", by.y="V1", all = F)
genexp2_norm_lnc_cons2$clu_all[is.na(genexp2_norm_lnc_cons2$clu_all)] = 0
genexp2_norm_lnc_all <- subset(genexp2_norm, str_detect(Genes, "align"))
genexp2_norm_lnc_all$clu_all[is.na(genexp2_norm_lnc_all$clu_all)] = 0
genexp2_norm_lnc_cluonly <- subset(genexp2_norm_lnc_all, clu_all != 0)

genexp2_lncRNAs <- subset(genexp2, str_detect(gene, "align"))

genexp2_lncRNAs_allichthyos <- merge(genexp2_lncRNAs, conservedLncs2, by.x="gene", by.y="V1", all = F)
genexp2_lncRNAs_allsphaeros <- merge(genexp2_lncRNAs, conservedLncs, by.x="gene", by.y="V1", all = F)
genexp2_lncRNAs_allsphaeros <- subset(genexp2_lncRNAs_allsphaeros, !(gene %in% conservedLncs2$V1))
genexp2_lncRNAs_onlySarc <- subset(genexp2_lncRNAs, !(gene %in% conservedLncs$V1))

### heatmaps of clustered lincRNAs
pal = (cividis(100))
superheat(as.matrix(genexp2_norm_lnc_cluonly[,2:21]),heat.pal=pal, membership.rows = genexp2_norm_lnc_cluonly$clu_all)

#boxplot of expression levels and cv of expression
boxplot(genexp2_lncRNAs_onlySarc$meanexp,genexp2_lncRNAs_allsphaeros$meanexp,genexp2_lncRNAs_allichthyos$meanexp, 
        outline = F, xaxt = "n", ylab = "mean expression (tpm)")
axis (1,at = 1:3, labels <- c("S.arctica-specific","conserved in \n Sphaeroforma spp", "conserved in other \n Ichthyosporeans"), cex.axis = 0.7)

boxplot(genexp2_lncRNAs_onlySarc$cv,genexp2_lncRNAs_allsphaeros$cv,genexp2_lncRNAs_allichthyos$cv, 
        outline = F, xaxt = "n", ylab = "coeffieicnt of variance of expression")
axis (1,at = 1:3, labels <- c("S.arctica-specific","conserved in \n Sphaeroforma spp", "conserved in other \n Ichthyosporeans"), cex.axis = 0.7)

###fraction of lnc expression per timepoint

totallncRNA <- colSums(genexp2[str_detect(genexp2$gene, "align"),2:21])/1e4

plot(totallncRNA[1:10],xlab = "time (h)",ylab = "total fraction of lncRNA expression (%)", type ="o", xaxt = "n", ylim = c(0,4))
axis(1,at = 1:10, labels <- 6*c(2:11))
lines(totallncRNA[11:20], lty = 4)
points(totallncRNA[11:20], pch = 16)
legend(7,4, lty = c(1,4), pch =c(1,16), legend = c("replicate 1","replicate 2"),bty = "n")


# plotting individual gene families ---------------------------------------

#Formins, Arp2/3/4 (actin nucleators)
genexp2_nucleators <- rbind(subset (genexp2_filt, str_detect(pfam,"FH2")), #formins,
                            subset (genexp2_filt, str_detect(gene,"Sarc4_g15431T")),#Arp2
                            subset (genexp2_filt, str_detect(gene,"Sarc4_g29717T")),#Arp3
                            subset (genexp2_filt, str_detect(gene,"Sarc4_g33563T")))#Arp4


pall = c(brewer.pal(6,"Blues"),brewer.pal(3,"Oranges"))
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))

plot(1:10,1:10,type = "n",ylim=c(0,50), xaxt = "n", xlab = "time (h)", ylab = "expression (tpm)", main = "ACTIN NUCLEATORS")
axis(1,at = 1:10, labels <- 6*c(2:11))
axis(4, at = c(0,10,20,30,40,50), labels = c(0,60,120,180,240,300))

for (i in c(1:6)) {
  lines(1:10,genexp2_nucleators[i,2:11], col = pall[i])
  points(1:10,genexp2_nucleators[i,2:11],pch=1, col = pall[i])
  lines(1:10,genexp2_nucleators[i,12:21],lty=4, col = pall[i])
  points(1:10,genexp2_nucleators[i,12:21],pch=16, col = pall[i])
} 
for (i in c(7:9)) {
  lines(1:10,1/6*genexp2_nucleators[i,2:11], col = pall[i])
  points(1:10,1/6*genexp2_nucleators[i,2:11],pch=1, col = pall[i])
  lines(1:10,1/6*genexp2_nucleators[i,12:21],lty=4, col = pall[i])
  points(1:10,1/6*genexp2_nucleators[i,12:21],pch=16, col = pall[i])
}     

legend(1,50, lty = 1, pch = 1, col = pall, bty="n", legend = c("Formin1",
                                                               "Formin2",
                                                               "Formin3",
                                                               "Formin4",
                                                               "Formin5",
                                                               "Formin6",
                                                               "Arp2*",
                                                               "Arp3*",
                                                               "Arp4*"))


#Septins,Cofilin,Profilin, MyoII,V (actin binding proteins)
genexp2_binding <- rbind(subset (genexp2_filt, str_detect(pfam,"Septin")), #septins,
                         subset (genexp2_filt, str_detect(gene,"Sarc4_g23768T")),#myoII
                         subset (genexp2_filt, str_detect(gene,"Sarc4_g13387T")),#myoV
                         subset (genexp2_filt, str_detect(gene,"Sarc4_g5267T")),#profilin
                         subset (genexp2_filt, str_detect(gene,"Sarc4_g6950T")))#cofilin

pall <- c(brewer.pal(8,"Blues")[c(5:8)],brewer.pal(3,"Oranges")[c(2:3)],brewer.pal(3,"Greens")[3],brewer.pal(3,"Reds")[3])

plot(1:10,1:10,type = "n",ylim=c(0,150), xaxt = "n", xlab = "time (h)", ylab = "expression (tpm)", main = "ACTIN-BINDING PROTEINS")
axis(1,at = 1:10, labels <- 6*c(2:11))
axis(4, at = c(0,40,80,120,160), labels = c(0,1000,2000,3000,4000))

for (i in c(1:6)) {
  lines(1:10,genexp2_binding[i,2:11], col = pall[i])
  points(1:10,genexp2_binding[i,2:11],pch=1, col = pall[i])
  lines(1:10,genexp2_binding[i,12:21],lty=4, col = pall[i])
  points(1:10,genexp2_binding[i,12:21],pch=16, col = pall[i])
} 
for (i in c(7:8)) {
  lines(1:10,1/25*genexp2_binding[i,2:11], col = pall[i])
  points(1:10,1/25*genexp2_binding[i,2:11],pch=1, col = pall[i])
  lines(1:10,1/25*genexp2_binding[i,12:21],lty=4, col = pall[i])
  points(1:10,1/25*genexp2_binding[i,12:21],pch=16, col = pall[i])
}     

legend(1,150, lty = 1, pch = 1, col = pall, bty="n", legend = c("Septin1",
                                                                "Septin2",
                                                                "Septin3",
                                                                "Septin4",
                                                                "Myosin II",
                                                                "Myosin V",
                                                                "Profilin / Chickadee *",
                                                                "Cofilin / Twinstar *"))

genexp2_adhesion <- rbind (subset (genexp2_filt, str_detect(gene,"Sarc4_g1310T")), #int_a
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g33805T")), #int_b 
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g2938T")), #pinch
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g25286T")), #parvin
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g20792T")), #paxilin
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g3705T")), #talin
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g11224T")), #vinculin
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g4850T")), #a-actinin
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g24897T")),#aardvark
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g24898T")),#aardvark
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g33617T")))#aardvark (not clustered)
pall = c(brewer.pal(8, "Spectral"),brewer.pal(3,"Greys"))


par(mfrow=c(1,1))
plot(1:10,1:10,type = "n",ylim=c(0,120), xaxt = "n", xlab = "time (h)", ylab = "expression (tpm)", main = "INTEGRIN ADHESOME AND AARDVARK")
axis(1,at = 1:10, labels <- 6*c(2:11))
axis(4, at = c(0,20,40,60,80,100,120), labels = c(0,100,200,300,400,500,600))
for (i in c(1:7,9,11)) {
  lines(1:10,genexp2_adhesion[i,2:11], col = pall[i])
  points(1:10,genexp2_adhesion[i,2:11],pch=1, col = pall[i])
  lines(1:10,genexp2_adhesion[i,12:21],lty=4, col = pall[i])
  points(1:10,genexp2_adhesion[i,12:21],pch=16, col = pall[i])
}  
for (i in c(8,10)) {
  lines(1:10,0.2*genexp2_adhesion[i,2:11], col = pall[i])
  points(1:10,0.2*genexp2_adhesion[i,2:11],pch=1, col = pall[i])
  lines(1:10,0.2*genexp2_adhesion[i,12:21],lty=4, col = pall[i])
  points(1:10,0.2*genexp2_adhesion[i,12:21],pch=16, col = pall[i])
}

legend(1,120, lty = 1, pch = 1, col = pall, bty="n", legend = c("integrin alpha",
                                                                "integrin beta",
                                                                "pinch",
                                                                "parvin",
                                                                "paxilin",
                                                                "talin",
                                                                "vinculin/alpha catenin",
                                                                "alpha actinin",
                                                                "aardvark1",
                                                                "aardvark2",
                                                                "aardvark3"))

###tubulins, 
genexp2_tubulins <- rbind (subset (genexp2_filt, str_detect(gene,"Sarc4_g31031T")), #tub a
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g29856T"))) #tub b 

pall = brewer.pal(3,"Greens")[c(2:3)]

par(mfrow=c(1,1))
plot(1:10,1:10,type = "n",ylim=c(0,1000), xaxt = "n", xlab = "time (h)", ylab = "mRNA abundance (tpm)", main = "tubulins")
axis(1,at = 1:10, labels <- 6*c(2:11))
for (i in c(1:2)) {
  lines(1:10,genexp2_tubulins[i,2:11], col = pall[i])
  points(1:10,genexp2_tubulins[i,2:11],pch=1, col = pall[i])
  lines(1:10,genexp2_tubulins[i,12:21],lty=4, col = pall[i])
  points(1:10,genexp2_tubulins[i,12:21],pch=16, col = pall[i])
}  

legend(1,1000, lty = 1, pch = 1, col = pall, bty="n", legend = c("alpha-tubulin",
                                                                 "beta-tubulin"))


genexp2_kinesins <- rbind (subset (genexp2_filt, str_detect(gene,"Sarc4_g8678T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g14701T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g23312T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g23711T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g23712T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g27337T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g29337T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g29515T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g31878T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g8472T")),
                           subset (genexp2_filt, str_detect(gene,"Sarc4_g8473T")))

pall = c(brewer.pal(9,"Oranges")[c(4:9)], brewer.pal(9,"Reds")[c(5:9)]) 

par(mfrow=c(1,1))
plot(1:10,1:10,type = "n",ylim=c(0,100), xaxt = "n", xlab = "time (h)", ylab = "mRNA abundance (tpm)", main = "kinesins")
axis(1,at = 1:10, labels <- 6*c(2:11))
for (i in c(1:11)) {
  lines(1:10,genexp2_kinesins[i,2:11], col = pall[i])
  points(1:10,genexp2_kinesins[i,2:11],pch=1, col = pall[i])
  lines(1:10,genexp2_kinesins[i,12:21],lty=4, col = pall[i])
  points(1:10,genexp2_kinesins[i,12:21],pch=16, col = pall[i])
}  

legend(1,100, lty = 1, pch = 1, col = pall, bty="n", legend = c("Kinesin 1",
                                                                "Kinesin 2",
                                                                "Kinesin 3",
                                                                "Kinesin 4",
                                                                "Kinesin 5",
                                                                "Kinesin 6",
                                                                "Kinesin 7",
                                                                "Kinesin 8",
                                                                "Kinesin 9",
                                                                "Kinesin 10",
                                                                "Kinesin 11"))

genexp2_ras <- rbind (subset (genexp2_filt, str_detect(gene,"Sarc4_g28650T")),
                      subset (genexp2_filt, str_detect(gene,"Sarc4_g29101T")),
                      subset (genexp2_filt, str_detect(gene,"Sarc4_g32829T")),
                      subset (genexp2_filt, str_detect(gene,"Sarc4_g5839T")))

pall = brewer.pal(4, "Blues")
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(1:10,1:10,type = "n",ylim=c(0,300), xaxt = "n", xlab = "time (h)", ylab = "expression (tpm)", main = "Rho GTPases")
axis(1,at = 1:10, labels <- 6*c(2:11))
axis(4, at = c(0,50,100,150,200,250,300), labels = c(0,200,400,600,800,1000,1200))
for (i in c(1,2,4)) {
  lines(1:10,genexp2_ras[i,2:11], col = pall[i])
  points(1:10,genexp2_ras[i,2:11],pch=1, col = pall[i])
  lines(1:10,genexp2_ras[i,12:21],lty=4, col = pall[i])
  points(1:10,genexp2_ras[i,12:21],pch=16, col = pall[i])
} 
lines(1:10,0.25*genexp2_ras[3,2:11], col = pall[3])
points(1:10,0.25*genexp2_ras[3,2:11],pch=1, col = pall[3])
lines(1:10,0.25*genexp2_ras[3,12:21],lty=4, col = pall[3])
points(1:10,0.25*genexp2_ras[3,12:21],pch=16, col = pall[3])

legend(1,300, lty = 1, pch = 1, col = pall, bty="n", legend = c("Rho1","Rho2","Rho3 *","Rho4"))


