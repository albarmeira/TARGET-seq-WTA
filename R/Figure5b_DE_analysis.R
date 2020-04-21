#High throughput scRNA-seq of primary myelofibrosis patient samples - 3'-TARGETseq dataset
#Assessment of transcriptional effects of different genetic subclones in the CD34+ compartment
#Differentially expressed genes in between each genetic subclone
#Author:Alba Rodriguez-Meira
#TARGET-seq : Unravelling intratumoral heterogeneity through high-sensitivity single-cell mutational analysis and parallel RNA-sequencing
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

###################################
#Perform differential expression analysis
################################### 

#JAK2_only vs all

#Select cells with JAK2 mutations but not other co-operating mutations
groupA<-as.character(cell.anno$cell.id[grep("JAK2_only",cell.anno$genotype_class1)]) 
groupA.m<-as.matrix(cell.count.use.log2[,colnames(cell.count.use.log2) %in% groupA])

#Select cells with JAK2 mutations and other co-operating mutations
groupB<-as.character(cell.anno$cell.id[grep("JAK2_only",cell.anno$genotype_class1,invert = TRUE)])
groupB.m<-as.matrix(cell.count.use.log2[,colnames(cell.count.use.log2) %in% groupB])
length(groupB)

###############################################################################
fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)

result.f<-data.frame()

for(i in 1:nrow(groupA.m)){
  
  n<-groupA.m[i,]
  m<-groupB.m[i,]
  n.f<-n[n>0]
  m.f<-m[m>0]
  
  x<-length(n)
  y<-length(m)
  
  x.1<-length(n.f)
  y.1<-length(m.f)
  x.2<-x-x.1
  y.2<-y-y.1
  
  m.m <- matrix(c(x.1, x.2, y.1, y.2), ncol = 2)
  fisher.test.p <-fisher.test(m.m)$p.value
  chi.test.p <-chisq.test(m.m)$p.value
  
  wilcox.p<-wilcox.test(n,m,paired = FALSE)$p.value
  
  all.p<-c(fisher.test.p,wilcox.p)
  fisher.p<-fishersMethod (all.p)
  
  my.meanA<-mean(groupA.m[i,])
  my.meanB<-mean(groupB.m[i,])
  my.log2.fc<-(my.meanB+1)-(my.meanA+1)
  my.test.f<-data.frame(gene=rownames(groupA.m)[i],meanA=my.meanA,meanB=my.meanB,
                        log2fc=my.log2.fc,nCellA=x,nCellB=y,expCellA=x.1,expCellB=y.1,
                        expFractionCellA=x.1/x,expFractionCellB=y.1/y,chisq.test.p=chi.test.p,fisher.test.p=fisher.test.p,wilcox=wilcox.p,fisher=fisher.p)
  result.f<-rbind(result.f,my.test.f)
}

#############################################################################
result.f<-result.f[order(result.f$fisher),] 
result.f$p.adjust<-p.adjust(result.f$fisher, method ="BH") # Calculate adjusted p-values
#############################################################################
DE_JAK2onlyvsall<-subset(result.f,p.adjust<0.1) #Select genes with adjusted p-values lower than 0.1
DE_JAK2onlyvsall<-subset(DE_JAK2onlyvsall,abs(log2fc)>1) #Select genes with log2 fold-change higher or lower than 1


