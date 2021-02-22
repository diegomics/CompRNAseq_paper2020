# Transcriptional responses are oriented towards different components of the rearing environment in two Drosophila sibling species.
# De Panis et al. 2020

# Inter species DE analysis

# Download/Install/Load required packages ##################################################

#install.packages("BiocManager")

#BiocManager::version()
#BiocManager::install("GenomicFeatures", version = "3.8")
#BiocManager::install("NOISeq", version = "3.8")
#BiocManager::install("DESeq2", version = "3.8")
#BiocManager::install("RColorBrewer", version = "3.8")
#BiocManager::install("ggplot2", version = "3.8")
#BiocManager::install("ggfortify", version = "3.8")
#BiocManager::install("rgl", version = "3.8")
#BiocManager::install("smacof", version = "3.8")
install.packages("export")

BiocManager::install("survival", version = "3.8", lib="/home/diego/R/x86_64-pc-linux-gnu-library/3.5")

library( GenomicFeatures )
library( NOISeq )
library( DESeq2 )
library( RColorBrewer )
library( ggplot2 )
library( ggfortify )
library( rgl )
library( smacof )
library( export )

setwd( "" )

# Import STAR/StringTie raw counts ##################################################

mycounts_ST <- as.matrix( read.csv( "gene_count_matrix.csv",
                                    header= TRUE, sep= ",", row.names= "gene_id" ) ) # raw count table by StringTie using mean read length  = 92
# mean read length was obtained using this bash script:
# awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sq/n-m*m);}'  *.fastq

head( mycounts_ST ) #rows (genes) and columns (treatments) are ordered by StringTie
dim( mycounts_ST ) #14680 genes

sample_info <- read.csv( "sample_info.csv",
                               header= TRUE, sep= "," )
head( sample_info )


# Get gene lenghts
# Import the GTF file used as input for STAR/StringTie
#txdb <- makeTxDbFromGFF("/home/diego/Dropbox/Projects/CompFly/RNASeq/3_Dmoj/StringTie/Op_A/Op_A.gtf",format="gtf")

# Collect the exons per gene id
#exons.list.per.gene <- exonsBy(txdb,by="gene")

# For each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
#exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

# Export lengths
#gene.lengths <- as.data.frame(do.call(rbind, exonic.gene.sizes))
#colnames(gene.lengths) <- NULL
#write.table(gene.lengths, "/home/diego/Dropbox/Projects/CompFly/RNASeq/3_Dmoj/StringTie/lengths", quote = F)

# Ordenar (en bash) la tabla de longitudes para que quede en el mismo orden que la tabla de conteos
#for line in $(cat gene_id_names); do grep -Fw $line lengths >> testc; done
gene_lengths <- read.csv( "/home/diego/Dropbox/Projects/CompFly/RNASeq/3_Dmoj/StringTie/gene_id.lengths",
                          header= TRUE, sep= " ", row.names= "gene_id" )

head(gene_lengths)
dim(gene_lengths)


# Check library sizes 
librarySizes <- colSums( mycounts_ST )
barplot( librarySizes, names= names( librarySizes ), las= 2, main= "Barplot of library sizes" )
# the low-coverage (2-rep treatments) have similar size but smaller than the rest of the tratments
# spliting the set allows better filtration by zeros and more accurate normalization

# normalization for plots only
coldata <- data.frame( Sample= sample_info$Sample,
                             Treatment= sample_info$Treatment.Species )

dds <- DESeqDataSetFromMatrix(countData = mycounts_ST,
                                    colData = coldata,
                                    design = ~ Treatment)
dds

rld.f <- rlog(dds, blind=FALSE)
head(assay(rld.f))

vst.f <- vst(dds, blind=FALSE)
head(assay(vst.f))

#
myfactors <- data.frame( treatments= c( "Op_B", "Op_B", "Op_K", "Op_K",
                                        "Op.Pr.Al_B", "Op.Pr.Al_B", "Op.Pr.Al_B", "Op.Pr.Al_K", "Op.Pr.Al_K", "Op.Pr.Al_K",
                                        "Op.Pr.Na_B", "Op.Pr.Na_B", "Op.Pr.Na_B", "Op.Pr.Na_K", "Op.Pr.Na_K", "Op.Pr.Na_K",
                                        "Tr.Pr.Al_B", "Tr.Pr.Al_B", "Tr.Pr.Al_B", "Tr.Pr.Al_K", "Tr.Pr.Al_K", "Tr.Pr.Al_K",
                                        "Tr.Pr.Na_B", "Tr.Pr.Na_B", "Tr.Pr.Na_B", "Tr.Pr.Na_K", "Tr.Pr.Na_K", "Tr.Pr.Na_K",
                                        "Tr_B", "Tr_B", "Tr_K", "Tr_K" ) )

mydata <- readData(data = mycounts_ST, factors = myfactors, length= gene_lengths$lengths)
#dat(mydata.biREP, type = "cd", norm = FALSE, refColumn = 1) #normalization needed
myTMM <- tmm(assayData(mydata)$exprs, long= gene_lengths$lengths, lc= 1, k=0.5)



# Subsets ##################################################

## Low nutrition subset = bireplicate (biREP) #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
mycounts.biREP <- as.matrix( read.csv( "/home/diego/Dropbox/Projects/CompFly/RNASeq/3_Dmoj/StringTie/gene_count_matrix_l92.csv",
                                          header= TRUE, sep= ",", row.names= "gene_id" )[, c( "Op_A", "Op_B",
                                                                                              "Op_D", "Op_E",
                                                                                              "Tr_A", "Tr_B",
                                                                                              "Tr_D", "Tr_E" )] )
head( mycounts.biREP )
dim( mycounts.biREP )
#
sample_info.biREP <- read.csv( "/home/diego/Dropbox/Projects/CompFly/RNASeq/3_Dmoj/DE/sample_info_biREP.csv",
                               header= TRUE, sep= "," )
head( sample_info.biREP )


# Check library sizes
librarySizes.biREP <- colSums( mycounts.biREP )
barplot( librarySizes.biREP, names= names( librarySizes.biREP ), las= 2, main= "Barplot of library sizes" )


## Native & 2X alkaloids subset = trireplicate (triREP) #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
mycounts.triREP <- as.matrix( read.csv( "/home/diego/Dropbox/Projects/CompFly/RNASeq/3_Dmoj/StringTie/gene_count_matrix_l92.csv",
                                         header= TRUE, sep= ",", row.names= "gene_id" )[, c( "Op.Pr.Al_A", "Op.Pr.Al_B", "Op.Pr.Al_C",
                                                                                             "Op.Pr.Al_D", "Op.Pr.Al_E", "Op.Pr.Al_F",
                                                                                             "Op.Pr.Na_A", "Op.Pr.Na_B", "Op.Pr.Na_C",
                                                                                             "Op.Pr.Na_D", "Op.Pr.Na_E", "Op.Pr.Na_F",
                                                                                             "Tr.Pr.Al_A", "Tr.Pr.Al_B", "Tr.Pr.Al_C",
                                                                                             "Tr.Pr.Al_D", "Tr.Pr.Al_E", "Tr.Pr.Al_F",
                                                                                             "Tr.Pr.Na_A", "Tr.Pr.Na_B", "Tr.Pr.Na_C",
                                                                                             "Tr.Pr.Na_D", "Tr.Pr.Na_E", "Tr.Pr.Na_F" )] )
head( mycounts.triREP )
dim( mycounts.triREP )

sample_info.triREP <- read.csv( "/home/diego/Dropbox/Projects/CompFly/RNASeq/3_Dmoj/DE/sample_info_triREP.csv",
                                header= TRUE, sep= "," )
head( sample_info.triREP )


# Check library sizes
librarySizes.triREP <- colSums( mycounts.triREP )
barplot( librarySizes.triREP, names= names( librarySizes.triREP ), las= 2, main= "Barplot of library sizes" )


# Prepare data for NOISeq ##################################################

## Low nutrition: biREP subset #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
myfactors.biREP <- data.frame( treatments= c( "Op_B", "Op_B", "Op_K", "Op_K",
                                              "Tr_B", "Tr_B", "Tr_K", "Tr_K" ) )

mydata.biREP <- readData(data = mycounts.biREP, factors = myfactors.biREP, length= gene_lengths$lengths)
#dat(mydata.biREP, type = "cd", norm = FALSE, refColumn = 1) #normalization needed
myTMM.biREP <- tmm(assayData(mydata.biREP)$exprs, long= gene_lengths$lengths, lc= 1, k=0.5)

## RLOG2/VST normalization (DESeq's) for plots only
# (performs a log2 scale transformation in a way that compensates for differences between samples for genes with low read count and also normalizes between samples for library size)
coldata_biREP <- data.frame( Sample= sample_info.biREP$Sample,
                     Treatment= sample_info.biREP$Treatment.Species )

dds_biREP <- DESeqDataSetFromMatrix(countData = mycounts.biREP,
                                    colData = coldata_biREP,
                                    design = ~ Treatment)
dds_biREP

rld_biREP.f <- rlog(dds_biREP, blind=FALSE)
head(assay(rld_biREP.f))

vst_biREP.f <- vst(dds_biREP, blind=FALSE)
head(assay(vst_biREP.f))


# Check distributions of samples using boxplots
boxplot( myTMM.biREP, xlab= "", ylab= "TMM(Counts)", las= 2 )
# Add a blue horizontal line that corresponds to the median logCPM
abline( h= median( as.matrix( myTMM.biREP ) ), col= "blue" )

# Exploratory PCA of the biREP subset
pcDat.biREP <- prcomp( t( myTMM.biREP ) )
# plot PCA
autoplot( pcDat.biREP, data= sample_info.biREP,
          colour= "Treatment", shape= "Species", size= 5 ) # + scale_shape_manual(values=c(21, 24)) + guides(fill = guide_legend(override.aes=list(shape=22)))
## Todos separan por especie.


## Native & 2X alkaloids: triREP subset #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
myfactors.triREP <- data.frame( treatments= c( "Op.Pr.Al_B", "Op.Pr.Al_B", "Op.Pr.Al_B", "Op.Pr.Al_K", "Op.Pr.Al_K", "Op.Pr.Al_K",
                                             "Op.Pr.Na_B", "Op.Pr.Na_B", "Op.Pr.Na_B", "Op.Pr.Na_K", "Op.Pr.Na_K", "Op.Pr.Na_K",
                                             "Tr.Pr.Al_B", "Tr.Pr.Al_B", "Tr.Pr.Al_B", "Tr.Pr.Al_K", "Tr.Pr.Al_K", "Tr.Pr.Al_K",
                                             "Tr.Pr.Na_B", "Tr.Pr.Na_B", "Tr.Pr.Na_B", "Tr.Pr.Na_K", "Tr.Pr.Na_K", "Tr.Pr.Na_K" ) )

mydata.triREP <- readData( data= mycounts.triREP, factors= myfactors.triREP, length= gene_lengths$lengths )
#dat(mydata.triREP, type = "cd", norm = FALSE, refColumn = 1) #normalization needed
myTMM.triREP <- tmm( assayData( mydata.triREP )$exprs, long= gene_lengths$lengths, lc= 1, k=0.5)


## RLOG2/VST normalization (DESeq's)
# (performs a log2 scale transformation in a way that compensates for differences between samples for genes with low read count and also normalizes between samples for library size)
coldata_triREP <- data.frame( Sample= sample_info.triREP$Sample,
                             Treatment= sample_info.triREP$Treatment.Species )

dds_triREP <- DESeqDataSetFromMatrix(countData = mycounts.triREP,
                                     colData = coldata_triREP,
                                     design = ~ Treatment)
dds_triREP

rld_triREP.f <- rlog(dds_triREP, blind=FALSE)
head(assay(rld_triREP.f))

vst_triREP.f <- vst(dds_triREP, blind=FALSE)
head(assay(vst_triREP.f))


# Check distributions of samples using boxplots
boxplot( myTMM.triREP, xlab= "", ylab= "TMM(Counts)", las= 2 )
# Add a blue horizontal line that corresponds to the median logCPM
abline( h= median( as.matrix( myTMM.triREP ) ), col= "blue" )

# Exploratory PCA of the triREP subset
pcDat.triREP <- prcomp( t( assay(rld_triREP.f) ) )
# plot PCA
autoplot( pcDat.triREP, data= sample_info.triREP, colour= "Treatment", shape= "Species", size= 5 )# + scale_shape_manual(values=c(21, 24)) + guides(fill = guide_legend(override.aes=list(shape=22)))
# solo la normalización de DESeq separa bien


# DE Analysis ##################################################

# (biREP) Op: Dbuz vs Dkoe #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

op_b.VS.op_k <- noiseqbio( mydata.biREP, k= 0.5, norm= "tmm",
                              factor= "treatments", conditions= c( "Op_B", "Op_K" ), lc= 1,
                              plot= FALSE, filter= 1 )

head( degenes( op_b.VS.op_k, q= 0.95, M= "down" ) )
# 0.00(NULL:11412) - 0.95(NULL:1885 | "up":1182 | "down":703) - 0.99(NULL:277 | "up":178 | "down":99)                 
#DE.plot(op_b.VS.op_k, q= 0.99, graphic= "expr", log.scale= TRUE) # Transcripts falling on the straight line (90° bisecting line) are equally expressed in both pools. Transcripts above the line are up-regulated in response to drought stress, while those below the line are down-regulated


# (biREP) Tr: Dbuz vs Dkoe #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
tr_b.VS.tr_k <- noiseqbio( mydata.biREP, k= 0.5, norm= "tmm",
                              factor= "treatments", conditions= c( "Tr_B", "Tr_K" ), lc= 1,
                              plot= FALSE, filter= 1 )

head( degenes( tr_b.VS.tr_k, q= 0.95, M= "down" ) )                  
# 0.00(NULL:11403) - 0.95(NULL:1176 | "up":823 | "down":353) - 0.99(NULL:139 | "up":89 | "down":50)                 
#DE.plot(tr_b.VS.tr_k, q= 0.99, graphic= "expr", log.scale= TRUE)

###

# (triREP) Op.Pr.Na: Dbuz vs Dkoe #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
op.pr.na_b.VS.op.pr.na_k <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                          factor= "treatments", conditions= c( "Op.Pr.Na_B", "Op.Pr.Na_K" ), lc= 1,
                                          plot= FALSE, filter= 1 )

head( degenes( op.pr.na_b.VS.op.pr.na_k, q= 0.95, M= "down" ) )
# 0.00(NULL:12194) - 0.95(NULL:4875 | "up":2463 | "down":2412) - 0.99(NULL:1498 | "up":832 | "down":666)                 
#DE.plot(op.pr.na_b.VS.op.pr.na_k, q= 0.99, graphic= "expr", log.scale= TRUE)


# (triREP) Op.Pr.Al: Dbuz vs Dkoe #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
op.pr.al_b.VS.op.pr.al_k <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                          factor= "treatments", conditions= c( "Op.Pr.Al_B", "Op.Pr.Al_K" ), lc= 1,
                                          plot= FALSE, filter= 1 )

head( degenes( op.pr.al_b.VS.op.pr.al_k, q= 0.95, M= "down" ) )
# 0.00(NULL:12338) - 0.95(NULL:4279 | "up":2582 | "down":1697) - 0.99(NULL:1184 | "up":840 | "down":344)                 
#DE.plot(op.pr.al_b.VS.op.pr.al_k, q= 0.99, graphic= "expr", log.scale= TRUE)


# (triREP) Tr.Pr.Na: Dbuz vs Dkoe #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
tr.pr.na_b.VS.tr.pr.na_k <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                          factor= "treatments", conditions= c( "Tr.Pr.Na_B", "Tr.Pr.Na_K" ), lc= 1,
                                          plot= FALSE, filter= 1 )

head( degenes( tr.pr.na_b.VS.tr.pr.na_k, q= 0.95, M= "down" ) )
# 0.00(NULL:12380) - 0.95(NULL:5036 | "up":2666 | "down":2370) - 0.99(NULL:1411 | "up":800 | "down":611)                 
#DE.plot(tr.pr.na_b.VS.tr.pr.na_k, q= 0.99, graphic= "expr", log.scale= TRUE)


# (triREP) Tr.Pr.Al: Dbuz vs Dkoe #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
tr.pr.al_b.VS.tr.pr.al_k <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                          factor= "treatments", conditions= c( "Tr.Pr.Al_B", "Tr.Pr.Al_K" ), lc= 1,
                                          plot= FALSE, filter= 1 )

head( degenes( tr.pr.al_b.VS.tr.pr.al_k, q= 0.95, M= "down" ) )  
# 0.00(NULL:12392) - 0.95(NULL:5863 | "up":3042 | "down":2821) - 0.99(NULL:1566 | "up":894 | "down":672)                 
#DE.plot(tr.pr.al_b.VS.tr.pr.al_k, q= 0.99, graphic= "expr", log.scale= TRUE)


# raw DE genes by subset ##################################################

# DE biREP (TMM normalized) #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
DE99_biREP.TMM <- subset( myTMM.biREP,
                          ifelse( row.names( myTMM.biREP ) %in% c(
                            rownames( degenes( op_b.VS.op_k, q= 0.99, M= NULL ) ),
                            rownames( degenes( tr_b.VS.tr_k, q= 0.99, M= NULL ) ) ),
                            yes= TRUE, no= FALSE ) )

dim( DE99_biREP.TMM ) #286 genes


# DE triREP (TMM normalized) #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
DE99_triREP.TMM <- subset( myTMM.triREP,
                           ifelse( row.names( myTMM.triREP ) %in% c(
                             rownames( degenes( op.pr.na_b.VS.op.pr.na_k, q= 0.99, M= NULL ) ),
                             rownames( degenes( op.pr.al_b.VS.op.pr.al_k, q= 0.99, M= NULL ) ),
                             rownames( degenes( tr.pr.na_b.VS.tr.pr.na_k, q= 0.99, M= NULL ) ),
                             rownames( degenes( tr.pr.al_b.VS.tr.pr.al_k, q= 0.99, M= NULL ) ) ),
                             yes= TRUE, no= FALSE ) )

dim( DE99_triREP.TMM ) #2952 genes

# ALL

DE99_ALL.TMM <- subset( myTMM,
                        ifelse( row.names( myTMM ) %in% c(
                          rownames( DE99_biREP.TMM ),
                          rownames( DE99_triREP.TMM ) ),
                          yes= TRUE, no= FALSE ) ) 

dim(DE99_ALL.TMM) #3000




# There are a few noisy DE genes...
# Analyze the set for a posteriori filter with the function "stats4filt":
stats4filt <- function(x){
  Min <- min(x, na.rm=TRUE)
  Max <- max(x, na.rm=TRUE)
  Median <- median(x, na.rm=TRUE)
  Mean <- mean(x, na.rm=TRUE)
  SD <- sd(x, na.rm=TRUE)
  
  return(c(Min=Min, Max=Max, Median=Median, Mean=Mean, SD=SD, rel1=SD/Mean, rel2=SD/Median, rel3=(SD*100)/Median, rel4=Max/Median ))
}

stats4filt_rslt.biREP.TMM <- as.data.frame(t(apply(DE99_biREP.TMM,1,stats4filt)))
stats4filt_rslt.triREP.TMM <- as.data.frame(t(apply(DE99_triREP.TMM,1,stats4filt)))

plot(stats4filt_rslt.triREP.TMM[,"Min"]) #br 60000
plot(stats4filt_rslt.biREP.TMM[,"Max"]) #br 2e5
plot(stats4filt_rslt.triREP.TMM[,"Median"]) #br 100000
plot(stats4filt_rslt.triREP.TMM[,"Mean"]) #br 125000
plot(stats4filt_rslt.triREP.TMM[,"SD"]) #br 80000 tr 3e5
plot(stats4filt_rslt.biREP.TMM[,"rel1"]) #br 2 ! tr 4 !!!
plot(stats4filt_rslt.biREP.TMM[,"rel2"]) #br 100 !!! tr 75 !!!
plot(stats4filt_rslt.biREP.TMM[,"rel3"]) #br 10000 !!! tr 1e5 !!!
plot(stats4filt_rslt.biREP.TMM[,"rel4"]) #br 500 !!! tr 5000 !!!


biREP.probNoise <- subset( DE99_biREP.TMM,
                           ifelse( row.names( DE99_biREP.TMM ) %in% c(
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel1"]>2),] ),
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel2"]>100),] ),
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel3"]>10000),] ),
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel4"]>500),] )
                           ),
                           yes = TRUE, no = FALSE))

dim(biREP.probNoise) #10

hmcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(288) # See "?brewer.pal" for more color schemes
H<-function(d) hclust(d, method="ward.D2")
heatmap( biREP.probNoise[, c(
  "Op_A","Op_B","Tr_A","Tr_B",
  "Op_D","Op_E","Tr_D","Tr_E" )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )


DE99_biREP.TMM.d <- subset( DE99_biREP.TMM,
                           ifelse( row.names( DE99_biREP.TMM ) %in% c(
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel1"]>2),] ),
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel2"]>100),] ),
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel3"]>10000),] ),
                             rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel4"]>500),] )
                           ),
                           yes = FALSE, no = TRUE))

dim(DE99_biREP.TMM.d) #276


#


triREP.probNoise <- subset( DE99_triREP.TMM,
                            ifelse( row.names( DE99_triREP.TMM ) %in% c(
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel1"]>4),] ),
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel2"]>75),] ),
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel3"]>1e5),] ),
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel4"]>5000),] )
                            ),
                            yes = TRUE, no = FALSE))

dim(triREP.probNoise) #25

hmcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(288) # See "?brewer.pal" for more color schemes
H<-function(d) hclust(d, method="ward.D2")
heatmap( triREP.probNoise[, c(
  "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C","Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
  "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C","Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
  "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F","Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
  "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F","Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )


DE99_triREP.TMM.d <- subset( DE99_triREP.TMM,
                            ifelse( row.names( DE99_triREP.TMM ) %in% c(
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel1"]>4),] ),
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel2"]>75),] ),
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel3"]>1e5),] ),
                              rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel4"]>5000),] )
                            ),
                            yes = FALSE, no = TRUE))

dim(DE99_triREP.TMM.d) #2927

#

DE99_ALL.TMM.d <- subset( myTMM,
                        ifelse( row.names( myTMM ) %in% c(
                          rownames( DE99_biREP.TMM.d ),
                          rownames( DE99_triREP.TMM.d ) ),
                          yes= TRUE, no= FALSE ) ) 

dim(DE99_ALL.TMM.d) #2972



# DE biREP - final dataset 

DE99_biREP.rld <- subset( assay(rld_biREP.f),
                          ifelse( row.names( assay(rld_biREP.f) ) %in% c(
                            rownames( DE99_biREP.TMM.d ) ),
                            yes= TRUE, no= FALSE ) )

DE99_biREP.vst <- subset( assay(vst_biREP.f),
                          ifelse( row.names( assay(vst_biREP.f) ) %in% c(
                            rownames( DE99_biREP.TMM.d ) ),
                            yes= TRUE, no= FALSE ) )

# exploratory PCA for DE biREP 
pcDat.DE99.biREP <- prcomp(t(DE99_biREP.vst))
# plot PCA
autoplot(pcDat.DE99.biREP, data = sample_info.biREP,
         colour="Treatment",
         shape="Species",
         size=5) +
  scale_color_manual(values = c("#ffc762","#e6f598"))

#scale_shape_manual(values=c(21, 24)) +
#  guides(fill = guide_legend(override.aes=list(shape=22)))

# exploratory heatmap for DE biREP 
hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" ) 
heatmap( DE99_biREP.vst[, c(
  "Op_A","Op_B","Tr_A","Tr_B",
  "Op_D","Op_E","Tr_D","Tr_E" )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )


# DE triREP - final dataset

DE99_triREP.rld <- subset( assay(rld_triREP.f),
                           ifelse( row.names( assay(rld_triREP.f) ) %in% c(
                             rownames( DE99_triREP.TMM.d ) ),
                             yes= TRUE, no= FALSE ) )

DE99_triREP.vst <- subset( assay(vst_triREP.f),
                           ifelse( row.names( assay(vst_triREP.f) ) %in% c(
                             rownames( DE99_triREP.TMM.d ) ),
                             yes= TRUE, no= FALSE ) )

# exploratory PCA for DE triREP 
pcDat.DE99.triREP <- prcomp(t(DE99_triREP.vst))
# plot PCA
autoplot(pcDat.DE99.triREP, data = sample_info.triREP,
         colour="Treatment",
         shape="Species",
         size=5) + 
  scale_color_manual(values = c("#9e0142","#f46d43","#3d85c6","#66c2a5"))
#scale_fill_manual(values = c("#FF1BB3","#A7FF5B","#99554D")) + scale_color_manual(values = c("black","white","orange"))


# exploratory heatmap for DE triREP  
hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
heatmap( DE99_triREP.vst[, c(
  "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C","Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
  "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C","Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
  "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F","Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
  "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F","Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

# ALL dataset - final

# ALL

DE99_ALL.rld <- subset( assay(rld.f),
                        ifelse( row.names( assay(rld.f) ) %in% c(
                          rownames( DE99_ALL.TMM.d ) ),
                          yes= TRUE, no= FALSE ) )            

DE99_ALL.vst <- subset( assay(vst.f),
                        ifelse( row.names( assay(vst.f) ) %in% c(
                          rownames( DE99_ALL.TMM.d ) ),
                          yes= TRUE, no= FALSE ) )            


hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" ) 
heatmap( DE99_ALL.vst[, c(
  "Op_A","Op_B","Tr_A","Tr_B",
  "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C","Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
  "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C","Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
  "Op_D","Op_E","Tr_D","Tr_E",
  "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F","Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
  "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F","Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F"  )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

pcDat.DE99.ALL <- prcomp(t(DE99_ALL.vst))
autoplot(pcDat.DE99.ALL, data = sample_info,
         colour="Treatment",
         shape="Species",
         size=5) 




# MDS ##################################################

d_ALL <- dist( t( DE99_ALL.vst ) )
MDS_ALL <- smacofSym(d_ALL, ndim=3) #
MDS_ALL$stress #

rgl.open()
rgl.bg(color = "white")

rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = "white" )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.8)
}

rgl.spheres( MDS_ALL$conf[1:3,1], MDS_ALL$conf[1:3,2], MDS_ALL$conf[1:3,3],
             color= "#9e0142", r= 0.02, lit= F ) #Op.Pr.Al_Dbuz
rgl.spheres( MDS_ALL$conf[7:9,1], MDS_ALL$conf[7:9,2], MDS_ALL$conf[7:9,3],
             color= "#f46d43", r= 0.02, lit= F ) #Op.Pr.Na_Dbuz

rgl.spheres( MDS_ALL$conf[4:6,1], MDS_ALL$conf[4:6,2], MDS_ALL$conf[4:6,3],
             color= "#9E0142", r= 0.02, lit= F ) #Op.Pr.Al_Dkoe
rgl.spheres( MDS_ALL$conf[10:12,1], MDS_ALL$conf[10:12,2], MDS_ALL$conf[10:12,3],
             color= "#F46D43", r= 0.02, lit= F ) #Op.Pr.Na_Dkoe

rgl.spheres( MDS_ALL$conf[13:14,1], MDS_ALL$conf[13:14,2], MDS_ALL$conf[13:14,3],
             color= "#ffc762", r= 0.02, lit= F ) #Op_Dbuz

rgl.spheres(  MDS_ALL$conf[15:16,1], MDS_ALL$conf[15:16,2], MDS_ALL$conf[15:16,3],
              color= "#ffc762", r= 0.02, lit= F ) #Op_Dkoe


rgl.spheres( MDS_ALL$conf[17:19,1], MDS_ALL$conf[17:19,2], MDS_ALL$conf[17:19,3],
             color= "#3d85c6", r= 0.02, lit= F ) #Tr.Pr.Al_Dbuz
rgl.spheres( MDS_ALL$conf[23:25,1], MDS_ALL$conf[23:25,2], MDS_ALL$conf[23:25,3],
             color= "#66c2a5", r= 0.02, lit= F ) #Tr.Pr.Na_Dbuz

rgl.spheres( MDS_ALL$conf[20:22,1], MDS_ALL$conf[20:22,2], MDS_ALL$conf[20:22,3],
             color= "#5E4FA2", r= 0.02, lit= F ) #Tr.Pr.Al_Dkoe
rgl.spheres( MDS_ALL$conf[26:28,1], MDS_ALL$conf[26:28,2], MDS_ALL$conf[26:28,3],
             color= "#66C2A5", r= 0.02, lit= F ) #Tr.Pr.Na_Dkoe

rgl.spheres( MDS_ALL$conf[29:30,1], MDS_ALL$conf[29:30,2], MDS_ALL$conf[29:30,3],
             color= "#e6f598", r= 0.02, lit= F ) #Tr_Dbuz

rgl.spheres(  MDS_ALL$conf[31:32,1], MDS_ALL$conf[31:32,2], MDS_ALL$conf[31:32,3],
              color= "#e6f598", r= 0.02, lit= F ) #Tr_Dkoe


ellipse_Dbuz <- ellipse3d( cov( cbind( MDS_ALL$conf[c( 1:3,7:9,13:14,17:19,23:25,29:30 ), 1],
                                             MDS_ALL$conf[c( 1:3,7:9,13:14,17:19,23:25,29:30 ), 2],
                                             MDS_ALL$conf[c( 1:3,7:9,13:14,17:19,23:25,29:30 ), 3]) ),
                                 centre= c(
                                   mean( MDS_ALL$conf[c( 1:3,7:9,13:14,17:19,23:25,29:30 ), 1]),
                                   mean( MDS_ALL$conf[c( 1:3,7:9,13:14,17:19,23:25,29:30 ), 2]),
                                   mean( MDS_ALL$conf[c( 1:3,7:9,13:14,17:19,23:25,29:30 ), 3]) ),
                                 scale= c( 1,1,1 ), level= 0.95 )

shade3d( ellipse_Dbuz, col= "#d9d9d9", alpha= 0.1, lit= F )
wire3d( ellipse_Dbuz, col= "#d9d9d9",  lit= F, alpha= 0.125 ) #los datos se concentran dentro de esa ellipse con una 

ellipse_Dkoe <- ellipse3d( cov( cbind( MDS_ALL$conf[c( 4:6,10:12,15:16,20:22,26:28,31:32 ), 1],
                                       MDS_ALL$conf[c( 4:6,10:12,15:16,20:22,26:28,31:32 ), 2],
                                       MDS_ALL$conf[c( 4:6,10:12,15:16,20:22,26:28,31:32 ), 3]) ),
                                 centre= c(
                                   mean( MDS_ALL$conf[c( 4:6,10:12,15:16,20:22,26:28,31:32 ), 1]),
                                   mean( MDS_ALL$conf[c( 4:6,10:12,15:16,20:22,26:28,31:32 ), 2]),
                                   mean( MDS_ALL$conf[c( 4:6,10:12,15:16,20:22,26:28,31:32 ), 3]) ),
                                 scale= c( 1,1,1 ), level= 0.95 )

shade3d( ellipse_Dkoe, col= "#666666", alpha= 0.1, lit= F )
wire3d( ellipse_Dkoe, col= "#666666",  lit= F, alpha= 0.125 ) #la ellipse representa una region con el 95% de confianza 


axes3d( edges= "bbox", labels= F, tick= F, box= F, expand= 1.025, color= "#252525" )
rgl.material(smooth = T, point_antialias = T, line_antialias = T)



rgl.postscript( "/home/diego/cubo_ALL.e.svg", fmt= "svg", drawText= T )


graph2svg(file = "/home/diego/Rchot", fun = NULL,
             aspectr = NULL, width = NULL, height = NULL, scaling = 100,
             bg = "white", colormodel = "rgb", cairo = TRUE)





# biREP #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

d_biREP <- dist( t( DE99_biREP.vst ) )

MDS_biREP <- smacofSym(d_biREP, ndim=3) #MDS ITERATIVE METRIC, NO PCA
MDS_biREP$stress #LOW STRESS, GOOD ANALYSIS!


rgl.open()
rgl.bg(color = "white")


rgl.spheres( MDS_biREP$conf[1:2,1], MDS_biREP$conf[1:2,2], MDS_biREP$conf[1:2,3],
             color= "#ffc762", r= 0.02, lit= F ) #Op_Dbuz

rgl.spheres(  MDS_biREP$conf[3:4,1], MDS_biREP$conf[3:4,2], MDS_biREP$conf[3:4,3],
             color= "#ffc762", r= 0.02, lit= F ) #Op_Dkoe


rgl.spheres( MDS_biREP$conf[5:6,1], MDS_biREP$conf[5:6,2], MDS_biREP$conf[5:6,3],
             color= "#e6f598", r= 0.02, lit= F ) #Tr_Dbuz

rgl.spheres(  MDS_biREP$conf[7:8,1], MDS_biREP$conf[7:8,2], MDS_biREP$conf[7:8,3],
             color= "#e6f598", r= 0.02, lit= F ) #Tr_Dkoe

# reg
x=c(MDS_biREP$conf[1:2,1],MDS_biREP$conf[5:6,1])
y=c(MDS_biREP$conf[1:2,3], MDS_biREP$conf[5:6,3])
z=c(MDS_biREP$conf[1:2,2], MDS_biREP$conf[5:6,2])

fit <- lm( y ~  x + z )
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
z.pred <- seq(min(z), max(z), length.out = grid.lines)
xz <- expand.grid( x = x.pred, z = z.pred)
y.pred <- matrix(predict(fit, newdata = xz), 
                 nrow = grid.lines, ncol = grid.lines)

rgl.surface(x.pred, z.pred, y.pred, color = "#d9d9d9", 
            alpha = 0.075, lit = FALSE)  

rgl.surface(x.pred, z.pred, y.pred, color = "#d9d9d9",
            alpha = 0.125, lit = FALSE, front = "lines", back = "lines")


x=c(MDS_biREP$conf[3:4,1],MDS_biREP$conf[7:8,1])
y=c(MDS_biREP$conf[3:4,3], MDS_biREP$conf[7:8,3])
z=c(MDS_biREP$conf[3:4,2], MDS_biREP$conf[7:8,2])

fit <- lm( y ~  x + z )
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
z.pred <- seq(min(z), max(z), length.out = grid.lines)
xz <- expand.grid( x = x.pred, z = z.pred)
y.pred <- matrix(predict(fit, newdata = xz), 
                 nrow = grid.lines, ncol = grid.lines)

rgl.surface(x.pred, z.pred, y.pred, color = "#666666", 
            alpha = 0.075, lit = FALSE)  

rgl.surface(x.pred, z.pred, y.pred, color = "#666666",
            alpha = 0.125, lit = FALSE, front = "lines", back = "lines")



axes3d( edges= "bbox", labels= F, tick= F, box= F, expand= 1.025, color= "#252525" )



# ellipses
ellipse_biREP_Dbuz <- ellipse3d( cov( cbind( MDS_biREP$conf[c( 1:2,5:6 ), 1],
                                             MDS_biREP$conf[c( 1:2,5:6 ), 2],
                                             MDS_biREP$conf[c( 1:2,5:6 ), 3]) ),
                                 centre= c(
                                   mean( MDS_biREP$conf[c( 1:2,5:6 ), 1]),
                                   mean( MDS_biREP$conf[c( 1:2,5:6 ), 2]),
                                   mean( MDS_biREP$conf[c( 1:2,5:6 ), 3]) ),
                                 scale= c( 1,1,1 ), level= 0.95 )

shade3d( ellipse_biREP_Dbuz, col= "#d9d9d9", alpha= 0.125, lit= F )

wire3d( ellipse_biREP_Dbuz, col= "#d9d9d9",  lit= F, alpha= 0.125 ) #

ellipse_biREP_Dkoe <- ellipse3d( cov( cbind( MDS_biREP$conf[c( 3:4,7:8 ), 1],
                                             MDS_biREP$conf[c( 3:4,7:8 ), 2],
                                             MDS_biREP$conf[c( 3:4,7:8 ), 3]) ),
                                 centre= c(
                                   mean( MDS_biREP$conf[c( 3:4,7:8 ), 1]),
                                   mean( MDS_biREP$conf[c( 3:4,7:8 ), 2]),
                                   mean( MDS_biREP$conf[c( 3:4,7:8 ), 3]) ),
                                 scale= c( 1,1,1 ), level= 0.95 )

shade3d( ellipse_biREP_Dkoe, col= "#666666", alpha= 0.0375, lit= F )
wire3d( ellipse_biREP_Dkoe, col= "#666666",  lit= F, alpha= 0.125 ) # ellipse 95% conf 


axes3d( edges= "bbox", labels= F, tick= F, box= F, expand= 1.025, color= "#252525" )



# triREP #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

d_triREP <- dist( t( DE99_triREP.vst ) )

MDS_triREP <- smacofSym( d_triREP, ndim= 3 ) #
MDS_triREP$stress #


rgl.open()
rgl.bg( color= "white" )

rgl.spheres( MDS_triREP$conf[1:3,1], MDS_triREP$conf[1:3,2], MDS_triREP$conf[1:3,3],
             color= "#9e0142", r= 0.02, lit= F ) #Op.Pr.Al_Dbuz
rgl.spheres( MDS_triREP$conf[7:9,1], MDS_triREP$conf[7:9,2], MDS_triREP$conf[7:9,3],
             color= "#f46d43", r= 0.02, lit= F ) #Op.Pr.Na_Dbuz

rgl.spheres( MDS_triREP$conf[4:6,1], MDS_triREP$conf[4:6,2], MDS_triREP$conf[4:6,3],
             color= "#9E0142", r= 0.02, lit= F ) #Op.Pr.Al_Dkoe
rgl.spheres( MDS_triREP$conf[10:12,1], MDS_triREP$conf[10:12,2], MDS_triREP$conf[10:12,3],
             color= "#F46D43", r= 0.02, lit= F ) #Op.Pr.Na_Dkoe


rgl.spheres( MDS_triREP$conf[13:15,1], MDS_triREP$conf[13:15,2], MDS_triREP$conf[13:15,3],
             color= "#3d85c6", r= 0.02, lit= F ) #Tr.Pr.Al_Dbuz
rgl.spheres( MDS_triREP$conf[19:21,1], MDS_triREP$conf[19:21,2], MDS_triREP$conf[19:21,3],
             color= "#66c2a5", r= 0.02, lit= F ) #Tr.Pr.Na_Dbuz

rgl.spheres( MDS_triREP$conf[16:18,1], MDS_triREP$conf[16:18,2], MDS_triREP$conf[16:18,3],
             color= "#5E4FA2", r= 0.02, lit= F ) #Tr.Pr.Al_Dkoe
rgl.spheres( MDS_triREP$conf[22:24,1], MDS_triREP$conf[22:24,2], MDS_triREP$conf[22:24,3],
             color= "#66C2A5", r= 0.02, lit= F ) #Tr.Pr.Na_Dkoe



x= c( MDS_triREP$conf[1:3,1], MDS_triREP$conf[7:9,1], MDS_triREP$conf[13:15,1], MDS_triREP$conf[19:21,1] )
y= c( MDS_triREP$conf[1:3,2], MDS_triREP$conf[7:9,2], MDS_triREP$conf[13:15,2], MDS_triREP$conf[19:21,2] )
z= c( MDS_triREP$conf[1:3,3], MDS_triREP$conf[7:9,3], MDS_triREP$conf[13:15,3], MDS_triREP$conf[19:21,3] )

fit <- lm( y ~  x + z )
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
z.pred <- seq(min(z), max(z), length.out = grid.lines)
xz <- expand.grid( x = x.pred, z = z.pred)
y.pred <- matrix(predict(fit, newdata = xz), 
                 nrow = grid.lines, ncol = grid.lines)

rgl.surface(x.pred, z.pred, y.pred, color = "#d9d9d9", 
            alpha = 0.075, lit = FALSE)  

rgl.surface(x.pred, z.pred, y.pred, color = "#d9d9d9",
            alpha = 0.125, lit = FALSE, front = "lines", back = "lines")


x= c( MDS_triREP$conf[4:6,1], MDS_triREP$conf[10:12,1], MDS_triREP$conf[16:18,1], MDS_triREP$conf[22:24,1] )
y= c( MDS_triREP$conf[4:6,2], MDS_triREP$conf[10:12,2], MDS_triREP$conf[16:18,2], MDS_triREP$conf[22:24,2] )
z= c( MDS_triREP$conf[4:6,3], MDS_triREP$conf[10:12,3], MDS_triREP$conf[16:18,3], MDS_triREP$conf[22:24,3] )

fit <- lm( y ~  x + z )
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
z.pred <- seq(min(z), max(z), length.out = grid.lines)
xz <- expand.grid( x = x.pred, z = z.pred)
y.pred <- matrix(predict(fit, newdata = xz), 
                 nrow = grid.lines, ncol = grid.lines)

rgl.surface(x.pred, z.pred, y.pred, color = "#666666", 
            alpha = 0.075, lit = FALSE)  

rgl.surface(x.pred, z.pred, y.pred, color = "#666666",
            alpha = 0.125, lit = FALSE, front = "lines", back = "lines")



axes3d( edges= "bbox", labels= F, tick= F, box= F, expand= 1.025, color= "#252525" )






ellipse_triREP_Dbuz <- ellipse3d( cov( cbind( MDS_triREP$conf[c( 1:3,7:9,13:15,19:21 ),1],
                                              MDS_triREP$conf[c( 1:3,7:9,13:15,19:21 ),2],
                                              MDS_triREP$conf[c( 1:3,7:9,13:15,19:21 ),3]) ),
                                  centre= c( mean( MDS_triREP$conf[c( 1:3,7:9,13:15,19:21 ),1]),
                                             mean( MDS_triREP$conf[c( 1:3,7:9,13:15,19:21 ),2]),
                                             mean( MDS_triREP$conf[c( 1:3,7:9,13:15,19:21 ),3]) ),
                                  scale= c( 1,1,1 ), level= 0.95 )

shade3d( ellipse_triREP_Dbuz, col= "#d9d9d9", alpha= 0.0375, lit= F )
wire3d( ellipse_triREP_Dbuz, col= "#d9d9d9",  lit= F, alpha= 0.125 ) # 


ellipse_triREP_Dkoe <- ellipse3d( cov( cbind( MDS_triREP$conf[c( 4:6,10:12,16:18,22:24 ),1],
                                              MDS_triREP$conf[c( 4:6,10:12,16:18,22:24 ),2],
                                              MDS_triREP$conf[c( 4:6,10:12,16:18,22:24 ),3]) ),
                                  centre= c( mean( MDS_triREP$conf[c( 4:6,10:12,16:18,22:24 ),1]),
                                             mean( MDS_triREP$conf[c( 4:6,10:12,16:18,22:24 ),2]),
                                             mean( MDS_triREP$conf[c( 4:6,10:12,16:18,22:24 ),3]) ),
                                  scale= c( 1,1,1 ), level= 0.95 )

shade3d( ellipse_triREP_Dkoe, col= "#666666", alpha= 0.0375, lit= F )
wire3d( ellipse_triREP_Dkoe, col= "#666666",  lit= F, alpha= 0.125 ) # 

axes3d( edges= "bbox", labels= F, tick= F, box= F, expand= 1.025, color= "#252525")


#TOP100

# We estimate the variance for each row in the logcounts matrix
countVar_ALL <- apply(DE99_ALL.vst, 1, var)
# Get the row numbers for the top 100 most variable genes
highVar_ALL <- order(countVar_ALL, decreasing=TRUE)[1:100]
# Subset logcounts matrix
top100_ALL <- DE99_ALL.vst[highVar_ALL,]

hm.TOP100 <- heatmap( top100_ALL[, c(
  "Op_A","Op_B","Tr_A","Tr_B",
  "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C","Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
  "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C","Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
  "Op_D","Op_E","Tr_D","Tr_E",
  "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F","Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
  "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F","Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F"  )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )









# biREP #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Enzimas
biREP.SLC <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/biREP_all.list_SLC", header = F, row.names = 1)
#
biREP.OXR <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/biREP_all.list_OXR", header = F, row.names = 1)
biREP.CE <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/biREP_all.list_CE", header = F, row.names = 1)
#
biREP.GST <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/biREP_all.list_GST", header = F, row.names = 1)
biREP.GT <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/biREP_all.list_GT", header = F, row.names = 1)
#
biREP.ABC <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/biREP_all.list_ABC", header = F, row.names = 1)


biREP.SLC.rlog <- subset( DE99_biREP.RLOG,
                          ifelse( row.names( DE99_biREP.RLOG ) %in% c( rownames( biREP.SLC ) ),
                                  yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" )
hm.biREP.SLC <- heatmap( biREP.SLC.rlog[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                             "Op_D","Op_E","Tr_D","Tr_E" )],
                         Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( biREP.SLC.rlog )[ hm.biREP.SLC$rowInd ] ),
            "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.SLC_rows",
            quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )
### keep order list
# hm.slc.db.rows <- rownames(SLC.DB_subset.rlog)[hm.slc.db$rowInd] #lista con el orden del dendrograma!!!!!
#https://sebastianraschka.com/Articles/heatmaps_in_r.html
#rownames(tess)[hm$rowInd] #lista con el orden del dendrograma!!!!!


biREP.OXR.rlog <- subset( DE99_biREP.RLOG,
                          ifelse( row.names( DE99_biREP.RLOG ) %in% c( rownames( biREP.OXR ) ),
                                  yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" )
hm.biREP.OXR <- heatmap( biREP.OXR.rlog[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                             "Op_D","Op_E","Tr_D","Tr_E" )],
                         Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( biREP.OXR.rlog )[ hm.biREP.OXR$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.OXR_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


biREP.CE.rlog <- subset( DE99_biREP.RLOG,
                         ifelse( row.names( DE99_biREP.RLOG ) %in% c( rownames( biREP.CE ) ),
                                 yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" )
hm.biREP.CE <- heatmap( biREP.CE.rlog[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                           "Op_D","Op_E","Tr_D","Tr_E" )],
                        Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( biREP.CE.rlog )[ hm.biREP.CE$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.CE_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


biREP.GST.rlog <- subset( DE99_biREP.RLOG,
                          ifelse( row.names( DE99_biREP.RLOG ) %in% c( rownames( biREP.GST ) ),
                                  yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" )
hm.biREP.GST <- heatmap( biREP.GST.rlog[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                             "Op_D","Op_E","Tr_D","Tr_E" )],
                         Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( biREP.GST.rlog )[ hm.biREP.GST$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.GST_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


biREP.GT.rlog <- subset( DE99_biREP.RLOG,
                         ifelse( row.names( DE99_biREP.RLOG ) %in% c( rownames( biREP.GT ) ),
                                 yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" )
hm.biREP.GT <- heatmap( biREP.GT.rlog[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                           "Op_D","Op_E","Tr_D","Tr_E" )],
                        Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( biREP.GT.rlog )[ hm.biREP.GT$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.GT_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


biREP.ABC.rlog <- subset( DE99_biREP.RLOG,
                          ifelse( row.names( DE99_biREP.RLOG ) %in% c( rownames( biREP.ABC ) ),
                                  yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" )
hm.biREP.ABC <- heatmap( biREP.ABC.rlog[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                             "Op_D","Op_E","Tr_D","Tr_E" )],
                         Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( biREP.ABC.rlog )[ hm.biREP.ABC$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.ABC_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


# updo

biREP.updo <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/biREP.updo.list", header = F, row.names = 1)

biREP.updo.rlog <- subset( DE99_biREP.RLOG,
                           ifelse( row.names( DE99_biREP.RLOG ) %in% c( rownames( biREP.updo ) ),
                                   yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" )
hm.biREP.updo <- heatmap( biREP.updo.rlog[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                               "Op_D","Op_E","Tr_D","Tr_E" )],
                          Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( biREP.updo.rlog )[ hm.biREP.updo$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.updo_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


# TOP100 

# We estimate the variance for each row in the logcounts matrix
countVar_biREP <- apply(DE99_biREP.vst, 1, var)
# Get the row numbers for the top 100 most variable genes
highVar_biREP <- order(countVar_biREP, decreasing=TRUE)[1:100]
# Subset logcounts matrix
hmDat_biREP <- DE99_biREP.vst[highVar_biREP,]

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" ) 
hm.biREP.100 <- heatmap( hmDat_biREP[, c( "Op_A","Op_B","Tr_A","Tr_B",
                                           "Op_D","Op_E","Tr_D","Tr_E" )],
                          Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( hmDat_biREP )[ hm.biREP.100$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.biREP.100_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


# HEATMAPS triREP #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#Enzimas
triREP.SLC <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/triREP_all.list_SLC", header = F, row.names = 1)
#
triREP.OXR <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/triREP_all.list_OXR", header = F, row.names = 1)
triREP.CE <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/triREP_all.list_CE", header = F, row.names = 1)
#
triREP.GST <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/triREP_all.list_GST", header = F, row.names = 1)
triREP.GT <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/triREP_all.list_GT", header = F, row.names = 1)
#
triREP.ABC <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/triREP_all.list_ABC", header = F, row.names = 1)


triREP.SLC.rlog <- subset( DE99_triREP.RLOG_den,
                           ifelse( row.names( DE99_triREP.RLOG_den ) %in% c( rownames( triREP.SLC ) ),
                                   yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
hm.triREP.SLC <- heatmap( triREP.SLC.rlog [,c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                               "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                               "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                               "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                               "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                               "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                               "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                               "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
                          Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( triREP.SLC.rlog )[ hm.triREP.SLC$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.SLC_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


triREP.OXR.rlog <- subset( DE99_triREP.RLOG_den,
                           ifelse( row.names( DE99_triREP.RLOG_den ) %in% c( rownames( triREP.OXR ) ),
                                   yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
hm.triREP.OXR <- heatmap( triREP.OXR.rlog [,c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                               "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                               "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                               "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                               "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                               "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                               "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                               "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
                          Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( triREP.OXR.rlog )[ hm.triREP.OXR$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.OXR_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


triREP.CE.rlog <- subset( DE99_triREP.RLOG_den,
                          ifelse( row.names( DE99_triREP.RLOG_den ) %in% c( rownames( triREP.CE ) ),
                                  yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
hm.triREP.CE <- heatmap( triREP.CE.rlog [,c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                             "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                             "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                             "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                             "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                             "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                             "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                             "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
                         Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( triREP.CE.rlog )[ hm.triREP.CE$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.CE_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


triREP.GST.rlog <- subset( DE99_triREP.RLOG_den,
                           ifelse( row.names( DE99_triREP.RLOG_den ) %in% c( rownames( triREP.GST ) ),
                                  yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
hm.triREP.GST <- heatmap( triREP.GST.rlog [,c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                               "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                               "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                               "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                               "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                               "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                               "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                               "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
                          Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( triREP.GST.rlog )[ hm.triREP.GST$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.GST_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


triREP.GT.rlog <- subset( DE99_triREP.RLOG_den,
                          ifelse( row.names( DE99_triREP.RLOG_den ) %in% c( rownames( triREP.GT ) ),
                                  yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
hm.triREP.GT <- heatmap( triREP.GT.rlog [,c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                             "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                             "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                             "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                             "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                             "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                             "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                             "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
                         Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( triREP.GT.rlog )[ hm.triREP.GT$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.GT_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


triREP.ABC.rlog <- subset( DE99_triREP.RLOG_den,
                           ifelse( row.names( DE99_triREP.RLOG_den ) %in% c( rownames( triREP.ABC ) ),
                                   yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
hm.triREP.ABC <- heatmap( triREP.ABC.rlog [,c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                               "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                               "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                               "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                               "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                               "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                               "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                               "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
                          Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( triREP.ABC.rlog )[ hm.triREP.ABC$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.ABC_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


#updo
triREP.updo <- read.table(file = "/home/diego/Dropbox/Projects/CompFly/del300319/triREP.updo.list", header = F, row.names = 1)

triREP.updo.rlog <- subset( DE99_triREP.RLOG_den,
                            ifelse( row.names( DE99_triREP.RLOG_den ) %in% c( rownames( triREP.updo ) ),
                                   yes= TRUE, no= FALSE ) )

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" )
hm.triREP.updo <- heatmap( triREP.updo.rlog [,c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                                 "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                                 "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                                 "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                                 "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                                 "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                                 "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                                 "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F" )],
                           Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( triREP.updo.rlog )[ hm.triREP.updo$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.updo_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


#TOP100

# We estimate the variance for each row in the logcounts matrix
countVar_triREP <- apply(DE99_triREP.TMM.d, 1, var)
# Get the row numbers for the top 100 most variable genes
highVar_triREP <- order(countVar_triREP, decreasing=TRUE)[1:100]
# Subset logcounts matrix
hmDat_triREP <- DE99_triREP.TMM.d[highVar_triREP,]

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 ) # See "?brewer.pal" for more color schemes
H <- function( d ) hclust( d, method= "ward.D2" ) 
hm.triREP.100 <- heatmap( hmDat_triREP[, c( "Op.Pr.Na_A","Op.Pr.Na_B","Op.Pr.Na_C",
                                            "Tr.Pr.Na_A","Tr.Pr.Na_B","Tr.Pr.Na_C",
                                            "Op.Pr.Al_A","Op.Pr.Al_B","Op.Pr.Al_C",
                                            "Tr.Pr.Al_A","Tr.Pr.Al_B","Tr.Pr.Al_C",
                                            "Op.Pr.Na_D","Op.Pr.Na_E","Op.Pr.Na_F",
                                            "Tr.Pr.Na_D","Tr.Pr.Na_E","Tr.Pr.Na_F",
                                            "Op.Pr.Al_D","Op.Pr.Al_E","Op.Pr.Al_F",
                                            "Tr.Pr.Al_D","Tr.Pr.Al_E","Tr.Pr.Al_F"  )],
                          Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )



write.table( rev( rownames( hmDat_triREP )[ hm.triREP.100$rowInd ] ),
             "/home/diego/Dropbox/Projects/CompFly/del300319/hm.triREP.100_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )

DE99_triREP.TMM.d["FBgn0136974",]

# T-SNE ####


