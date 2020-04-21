#High throughput scRNA-seq of primary myelofibrosis patient samples - 3'-TARGETseq dataset
#Assessment of transcriptional effects of different genetic subclones in the CD34+ compartment
#Comparison all major subclones found in patient samples
#Author:Alba Rodriguez-Meira
#TARGET-seq : Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing
#PMID:30765193
###############################################################################

###################################
load(file="HT_TARGETseq_counts.qc.normalized.rdata")
load(file="HT_TARGETseq_colData.qc.rdata")
my.db<- gmtPathways("Human_signatures.gmt") #Human signatures - including cell cycle related genes
###################################

###################################
#Select cells from patients with major genetic subclones (over 5%)
###################################
cell.anno<-colData.qc

cell.anno<-subset(cell.anno,genotype.class=="realclone") #Remove undetermined cells or clones with less than 10 cells in the overall dataset (minor)

#Perform classification of genotypes

cell.anno$genotype_class1<-"NA"
JAK2_only<-c("\\bJAK2_HET\\b","\\bJAK2_HOM\\b")
cell.anno$genotype_class1[grep(paste(JAK2_only,collapse="|"),cell.anno$genotypesall)]<-"JAK2_only"

WT_normal<-c("\\bWT\\b")
cell.anno$genotype_class1[grep(paste(WT_normal,collapse="|"),cell.anno$genotypesall)]<-"WT_normal"

WT_patient<-c("\\bJAK2_WT\\b")
cell.anno$genotype_class1[grep(paste(WT_patient,collapse="|"),cell.anno$genotypesall)]<-"WT_patient"

JAK2_spliceosome<-c("JAK2_HET_SRSF2_HET","JAK2_HET_U2AF1_HET","JAK2_HOM_CBL_p404_HET_SRSF2_HET","JAK2_HOM_SF3B1_HET","JAK2_HOM_SRSF2_HET",
                    "JAK2_HOM_U2AF1_HET")
cell.anno$genotype_class1[grep(paste(JAK2_spliceosome,collapse="|"),cell.anno$genotypesall)]<-"JAK2_spliceosome"


JAK2_spliceosome_epigenetic<-c("JAK2_HET_U2AF1_HET_ASXL1_p910_HET","JAK2_HET_U2AF1_HET_TET2_pI1105_HET_ASXL1_p910_HET",
                               "JAK2_HET_U2AF1_HET_ASXL1_p897_HET")
cell.anno$genotype_class1[grep(paste(JAK2_spliceosome_epigenetic,collapse="|"),cell.anno$genotypesall)]<-"JAK2_spliceosome_epigenetic"

JAK2_epigenetic<-c("JAK2_HOM_ASXL1_p644_HET","JAK2_HOM_TET2_p1612_HET")
cell.anno$genotype_class1[grep(paste(JAK2_epigenetic,collapse="|"),cell.anno$genotypesall)]<-"JAK2_epigenetic"

cell.anno<-subset(cell.anno,genotype_class1!="NA") #2734 cells considered for this analysis
###################################

###################################
#Subset counts from selected cells 
cell.count<-as.data.frame(counts.qc.normalized)
cell.count<-cell.count[as.character(cell.anno$cell.id)]
###################################

###################################
#Remove ERCC from analysis
cell.count.no.ERCC <-cell.count[-c(grep("ERCC-", rownames(cell.count))),]
###################################

###################################
#Remove genes expressed in less than 5 cells
genes_with_expressing_cells<-5
cells.per.genes<-as.data.frame(apply(cell.count.no.ERCC>1,1,sum))
colnames(cells.per.genes)<-c("cell.exp.genes")
genes<-rownames(cells.per.genes[cells.per.genes$cell.exp.genes>5,,drop=FALSE])
cell.count.use<-as.matrix(cell.count.no.ERCC[genes,])
#############################################################

###################################
#Perform log2 transformation
cell.count.use.log2<-cell.count.use
min_expr<-1
cell.count.use.log2[cell.count.use.log2 < min_expr]<- 0
cell.count.use.log2[cell.count.use.log2 >= min_expr]<-log2(cell.count.use.log2[cell.count.use.log2 >= min_expr])
#############################################################

#############################################################
############Remove cell cycle effect#########################
#############################################################
###################################
#Extract cell cycle related genes 

g2m.genes<-as.character(my.db[["G2M_Tirosh_et_al_Science_2016"]])
s_phase.genes<-as.character(my.db[["S_phase_Tirosh_et_al_Science_2016"]])

my.g2m.exp.log2<-as.data.frame(cell.count.use.log2[rownames(cell.count.use.log2) %in% g2m.genes,])
my.g2m.exp.log2<-as.data.frame(t(my.g2m.exp.log2))
my.g2m.exp.log2$averageG2Mphase<-rowMeans(my.g2m.exp.log2)
my.g2m.exp.log2$Cell<-as.character(rownames((my.g2m.exp.log2)))

my.data.g2M<-merge(cell.anno,my.g2m.exp.log2,by.x="cell.id",by.y="Cell",all=T)

my.Sphase.exp.log2<-as.data.frame(cell.count.use.log2[rownames(cell.count.use.log2) %in% s_phase.genes,])
my.Sphase.exp.log2<-as.data.frame(t(my.Sphase.exp.log2))
my.Sphase.exp.log2$averageSphase<-rowMeans(my.Sphase.exp.log2)
my.Sphase.exp.log2$Cell<-as.character(rownames((my.Sphase.exp.log2)))

my.data.Sphase<-merge(cell.anno,my.Sphase.exp.log2,by.x="cell.id",by.y="Cell",all=T)

my.data.cell.cycle<-merge(my.data.Sphase,my.data.g2M,by.x="cell.id",by.y="cell.id",all=T)

###################################
#Remove unwanted sources of variation

FM<-cell.count.use.log2[,as.character(my.data.cell.cycle$cell.id)]

set.seed(1)
X.model_mat <- Matrix::sparse.model.matrix(as.formula("~averageSphase+averageG2Mphase"),
                                           data = my.data.cell.cycle, drop.unused.levels = TRUE)
fit <- limma::lmFit(FM, X.model_mat)
beta <- fit$coefficients[, -1, drop = FALSE]
beta[is.na(beta)] <- 0
FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
FM <- FM[!is.na(row.names(FM)), ]

cell.count.use.log2.nocellcycle<-as.data.frame(as.matrix(FM))
#########################################################

#########################################################
###########Remove batch effect###########################
#########################################################

my_genotype<-as.character(cell.anno$genotype_class1)
my.design<-model.matrix(~my_genotype)
rownames(my.design)<-colnames(cell.count.use.log2.nocellcycle)

batch.x<-as.character(cell.anno$batch)
batch.y<-as.character(cell.anno$donor)

counts.no.batch = removeBatchEffect(cell.count.use.log2.nocellcycle,batch=batch.x,batch2=batch.y,design=my.design)
####################################################################

#########################################################
#############Get variable genes##########################
#########################################################

my.mean<-rowMeans(cell.count.use[,c(1:ncol(cell.count.use))])
vars.S <- rowSds(cell.count.use)
my.cv <- vars.S / my.mean

###############################
#Model CV vs Mean  

n.m<-data.frame(mean=log2(my.mean),cv=log2(my.cv))
rownames(n.m)<-rownames(cell.count.use)
n.m<-na.omit(n.m)

lowess <- loess(n.m$cv~n.m$mean,span=0.65)

n.m$predicted.cv<-predict(lowess,n.m$mean)
n.m$usedGene<-F
n.m$usedGene[n.m$predicted.cv < n.m$cv & n.m$mean >= 0]<-T
n.m.used<-subset(n.m,usedGene ==T)

plot(n.m$mean,n.m$cv,pch=19,col="black",xlab="Mean of normalized reads, Log2",ylab="Coefficient of variation (CV), Log2",cex=0.2,
     main="")
xl <- seq(min(n.m$mean), max(n.m$mean), length.out=100)
points(n.m.used$mean,n.m.used$cv,col="cyan",pch=19,cex=0.5)
lines(xl, predict(lowess,xl), col='red', lwd=2)
text(5,4,label=paste(nrow(n.m.used)," genes",sep=""))

############################################
hv.genes<-as.character(rownames(n.m.used)) #3286 hv genes
final.matrix<-subset(counts.no.batch,rownames(counts.no.batch) %in% hv.genes)

############################################
##############PCA analysis##################
############################################

set.seed(1)
res_pca<-prcomp_irlba(t(final.matrix), n=50, center = TRUE,scale. = TRUE)

prop_varex <- res_pca$sdev^2/sum(res_pca$sdev^2)
p <- qplot(1:length(prop_varex), prop_varex, alpha = I(0.5)) + theme(legend.position = "top",
                                                                     legend.key.height = grid::unit(0.35, "in")) +
  xlab("components") + ylab("Standard deviation of variance\n for each component")
p

############################################
#Run TSNE

n.dims.use<-30
my.pca<-as.matrix(res_pca$x[,1:n.dims.use])
set.seed(seed = 100)
n.perplexity<-30
tsne_out <- Rtsne(my.pca,check_duplicates = FALSE,perplexity = n.perplexity)

my.results<-data.frame(Cell=colnames(cell.count.use.log2),tsne_out$Y)
my.dat<-merge(my.results,cell.anno,by.x="Cell",by.y="cell.id",all=T)

############################################
#Plot Figure

ggplot(data=my.dat,aes(x=X1,y=X2,colour=genotype_class1))+geom_point(size=1.5)+
  theme_bw()+
  geom_point(size=1.5,shape=21,color='black',stroke=0.05,alpha=0.25)+
  scale_colour_manual(name="legend", values=c("purple","blue","olivedrab1","green4","grey90","grey50"))+
  theme(text = element_text(size=10))+
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background = element_blank())+
  theme(legend.position = "none")

ggsave(filename="Figure5a.pdf",width=6,height = 6)
