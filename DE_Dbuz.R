# Transcriptional responses are oriented towards different components of the rearing environment in two Drosophila sibling species.
# De Panis et al. 2020

# Intra species differential gene expression analysis for D. buzzatii (the same script was used for D. koepferae)

# Download/Install/Load required packages ##################################################

#install.packages( "BiocManager" )

#BiocManager::version()
#BiocManager::install( "GenomicFeatures" )
#BiocManager::install( "NOISeq" )
#BiocManager::install( "DESeq2" )
#BiocManager::install( "RColorBrewer" )
#BiocManager::install( "ggplot2" )
#BiocManager::install( "ggfortify" )
#BiocManager::install( "rgl" )
#BiocManager::install( "smacof" )

library( GenomicFeatures )
library( NOISeq )
library( DESeq2 )
library( RColorBrewer )
library( ggplot2 )
library( ggfortify )
library( rgl )
library( smacof )


setwd( "WORKING_DIR" )


# Import STAR/StringTie raw counts ##################################################

mycounts_ST <- as.matrix( read.csv( "Db.gene_count_matrix.csv",
                                    header= TRUE, sep= ",", row.names= "gene_id" ) ) # raw count table by StringTie using mean read length  = 92

# mean read length was obtained using this bash script:
# awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sq/n-m*m);}'  *.fastq

head( mycounts_ST ) #rows (genes) and columns (treatments) are ordered by StringTie
dim( mycounts_ST ) #13567

sample_info <- read.csv( "Db.sample_info.csv",
                         header= TRUE, sep= "," )
head( sample_info )

# Import gene lenghts
gene_lengths <- read.csv( "Db.gene_id.lengths",
                          header= TRUE, sep= " ", row.names= "gene_id" )

head(gene_lengths)
dim(gene_lengths)

# Define factors and TMM normalization for later 
# Db=Drosophila buzzatii; Op=Opuntia sulphurea; Tr=Trichocereus terscheckii; Lo=Low nutrition; Na=Native; 2X=2X alkaloids
myfactors <- data.frame( treatments= c( "Op.Lo_Db", "Op.Lo_Db",
                                        "Op.2X_Db", "Op.2X_Db", "Op.2X_Db",
                                        "Op.Na_Db", "Op.Na_Db", "Op.Na_Db",
                                        "Tr.2X_Db", "Tr.2X_Db", "Tr.2X_Db",
                                        "Tr.Na_Db", "Tr.Na_Db", "Tr.Na_Db",
                                        "Tr.Lo_Db", "Tr.Lo_Db" ) )

mydata <- readData(data = mycounts_ST, factors = myfactors, length= gene_lengths$lengths)
myTMM <- tmm(assayData(mydata)$exprs, long= gene_lengths$lengths, lc= 1, k=0.5)


# Check library sizes 
librarySizes <- colSums( mycounts_ST )
barplot( librarySizes, names= names( librarySizes ), las= 2, main= "Barplot of library sizes" )
# the low-coverage (2-rep treatments) have similar size but smaller than the rest of the tratments
# spliting the set in two subsets allows better filter by zeros & more accurate normalization


# Subsets ##################################################

# Low nutrition & Native (genotypes A, B) subset = bi-replicate = biREP #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
mycounts.biREP <- as.matrix( read.csv( "Db.gene_count_matrix.csv",
                                          header= TRUE, sep= ",", row.names= "gene_id" )[, c( "Op.Na_A", "Op.Na_B",
                                                                                              "Tr.Na_A", "Tr.Na_B",
                                                                                              "Op.Lo_A", "Op.Lo_B",
                                                                                              "Tr.Lo_A", "Tr.Lo_B" )] )
head( mycounts.biREP )
dim( mycounts.biREP )
#
sample_info.biREP <- read.csv( "Db.sample_info_biREP.csv",
                               header= TRUE, sep= "," )
head( sample_info.biREP )


# Check again library sizes
librarySizes.biREP <- colSums( mycounts.biREP )
barplot( librarySizes.biREP, names= names( librarySizes.biREP ), las= 2, main= "Barplot of library sizes" )


# Native & 2X alkaloids (genotypes A, B, C) subset = tri-replicate = triREP #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
mycounts.triREP <- as.matrix( read.csv( "Db.gene_count_matrix.csv",
                                           header= TRUE, sep= ",", row.names= "gene_id" )[, c( "Op.2X_A", "Op.2X_B", "Op.2X_C",
                                                                                               "Op.Na_A", "Op.Na_B", "Op.Na_C",
                                                                                               "Tr.2X_A", "Tr.2X_B", "Tr.2X_C",
                                                                                               "Tr.Na_A", "Tr.Na_B", "Tr.Na_C" )] )
head( mycounts.triREP )
dim( mycounts.triREP )

sample_info.triREP <- read.csv( "Db.sample_info_triREP.csv",
                                header= TRUE, sep= "," )
head( sample_info.triREP )


# Check again library sizes
librarySizes.triREP <- colSums( mycounts.triREP )
barplot( librarySizes.triREP, names= names( librarySizes.triREP ), las= 2, main= "Barplot of library sizes" )


# Prepare data for NOISeq and explore data clustering ##################################################

# Cactus & nutrition: biREP subset

# Check distributions of samples using boxplots (before normalization)
boxplot( mycounts.biREP, xlab= "", ylab= "Counts", las= 2 )
# Add a blue horizontal line that corresponds to the median logCPM
abline( h= median( as.matrix( mycounts.biREP ) ), col= "blue" )

myfactors.biREP <- data.frame( treatments= c( "Op.Na", "Op.Na",
                                              "Tr.Na", "Tr.Na",
                                              "Op.Lo", "Op.Lo",
                                              "Tr.Lo", "Tr.Lo" ) )

mydata.biREP <- readData(data = mycounts.biREP, factors = myfactors.biREP, length= gene_lengths$lengths)
myTMM.biREP <- tmm(assayData(mydata.biREP)$exprs, long= gene_lengths$lengths, lc= 1, k=0.5)

# In addition use VST normalization (DESeq's) for plotting
# (performs a log2 scale transformation in a way that compensates for differences between samples for genes with low read count and also normalizes between samples for library size)
coldata_biREP <- data.frame( Sample= sample_info.biREP$Sample,
                             Treatment= sample_info.biREP$Treatment )
head(coldata_biREP)

dds_biREP <- DESeqDataSetFromMatrix(countData = mycounts.biREP,
                                    colData = coldata_biREP,
                                    design = ~ Treatment)
dds_biREP

vst_biREP.f <- vst(dds_biREP, blind=FALSE)
head(assay(vst_biREP.f))

# Check distributions of normalized samples using boxplots
boxplot( myTMM.biREP, xlab= "", ylab= "NormCounts", las= 2 )
boxplot( assay(vst_biREP.f), xlab= "", ylab= "NormCounts", las= 2 )

# Exploratory PCA of the biREP subset  
pcDat.biREP <- prcomp( t( myTMM.biREP ) )
# plot PCA
autoplot( pcDat.biREP, data= sample_info.biREP,
          colour= "Protein", shape= "Cactus", size= 5 )


# Cactus & alkaloids: triREP subset

# Check distributions of samples using boxplots
boxplot( mycounts.triREP, xlab= "", ylab= "Counts", las= 2 )
# Add a blue horizontal line that corresponds to the median logCPM
abline( h= median( as.matrix( myTMM.triREP ) ), col= "blue" )

myfactors.triREP <- data.frame( treatments= c( "Op.2X", "Op.2X", "Op.2X",
                                               "Op.Na", "Op.Na", "Op.Na",
                                               "Tr.2X", "Tr.2X", "Tr.2X",
                                               "Tr.Na", "Tr.Na", "Tr.Na" ) )

mydata.triREP <- readData( data= mycounts.triREP, factors= myfactors.triREP, length= gene_lengths$lengths )
myTMM.triREP <- tmm( assayData( mydata.triREP )$exprs, long= gene_lengths$lengths, lc= 1, k=0.5)

# In addition use VST normalization (DESeq's) for plotting
coldata_triREP <- data.frame( Sample= sample_info.triREP$Sample,
                             Treatment= sample_info.triREP$Treatment )
head(coldata_triREP)

dds_triREP <- DESeqDataSetFromMatrix(countData = mycounts.triREP,
                                    colData = coldata_triREP,
                                    design = ~ Treatment)
dds_triREP

vst_triREP.f <- vst(dds_triREP, blind=FALSE)
head(assay(vst_triREP.f))


# Check distributions of normalized samples using boxplots
boxplot( myTMM.triREP, xlab= "", ylab= "NormCounts", las= 2 )
boxplot( assay(vst_triREP.f), xlab= "", ylab= "NormCounts", las= 2 )

# Exploratory PCA of the triREP subset
pcDat.triREP <- prcomp( t( myTMM.triREP ) )
# plot PCA
autoplot( pcDat.triREP, data= sample_info.triREP, colour= "Alk.Conc.", shape= "Cactus", size= 5 )


# Differential expression analysis ##################################################

# (biREP) Low nutrition: O. sulphurea 'Low nutrition' vs T. terscheckii 'Low nutrition' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
Op.Lo.VS.Tr.Lo <- noiseqbio( mydata.biREP, k= 0.5, norm= "tmm",
                       factor= "treatments", conditions= c( "Op.Lo", "Tr.Lo" ), lc= 1,
                       plot= F, filter= 1 )

head( degenes( Op.Lo.VS.Tr.Lo, q= 0.99, M= "up" ) ) # "up" for upregulated, "down" for downregulated, NULL for all DE
#DE.plot(Op.Lo.VS.Tr.Lo, q= 0.99, graphic= "expr", log.scale= TRUE) # Transcripts falling on the straight line (90° bisecting line) are equally expressed in both pools. Transcripts above the line are up-regulated in response to drought stress, while those below the line are down-regulated

# (biREP) O. sulphurea: 'Low nutrition' vs 'Native' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
Op.Lo.VS.Op.Na <- noiseqbio( mydata.biREP, k= 0.5, norm= "tmm",
                             factor= "treatments", conditions= c( "Op.Lo", "Op.Na" ), lc= 1,
                             plot= F, filter= 1 )

head( degenes( Op.Lo.VS.Op.Na, q= 0.99, M= "down" ) )

# (biREP) T. terscheckii: 'Low nutrition' vs 'Native' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
Tr.Lo.VS.Tr.Na <- noiseqbio( mydata.biREP, k= 0.5, norm= "tmm",
                             factor= "treatments", conditions= c( "Tr.Lo", "Tr.Na" ), lc= 1,
                             plot= F, filter= 1 )

head( degenes( Tr.Lo.VS.Tr.Na, q= 0.99, M= "down" ) )


###

# (triREP) 'Native': O. sulhpurea vs T. terscheckii #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
Op.Na.VS.Tr.Na <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                   factor= "treatments", conditions= c( "Op.Na", "Tr.Na" ), lc= 1,
                                   plot= F, filter= 1 )

head( degenes( Op.Na.VS.Tr.Na, q= 0.99, M= "down" ) )

# (triREP) '2X alkaloids': O. sulhpurea vs T. terscheckii #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
Op.2X.VS.Tr.2X <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                   factor= "treatments", conditions= c( "Op.2X", "Tr.2X" ), lc= 1,
                                   plot= F, filter= 1 )

head( degenes( Op.2X.VS.Tr.2X, q= 0.99, M= "down" ) )

# (triREP) O. sulphurea: 'Native' vs '2X alkaloids' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
Op.Na.VS.Op.2X <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                   factor= "treatments", conditions= c( "Op.Na", "Op.2X" ), lc= 1,
                                   plot= F, filter= 1 )

head( degenes( Op.Na.VS.Op.2X, q= 0.99, M= "down" ) )

# (triREP) T. terscheckii: 'Native' vs '2X alkaloids' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
Tr.Na.VS.Tr.2X <- noiseqbio( mydata.triREP, k= 0.5, norm= "tmm",
                                   factor= "treatments", conditions= c( "Tr.Na", "Tr.2X" ), lc= 1,
                                   plot= FALSE, filter= 1 )

head( degenes( Tr.Na.VS.Tr.2X, q= 0.99, M= "down" ) )

# raw DE genes by subset ##################################################

# DE biREP (TMM normalized)

DE99_biREP.TMM <- subset( myTMM.biREP,
                           ifelse( row.names( myTMM.biREP ) %in% c(
                             rownames( degenes( Op.Lo.VS.Tr.Lo, q= 0.99, M= NULL ) ),
                             rownames( degenes( Op.Lo.VS.Op.Na, q= 0.99, M= NULL ) ),
                             rownames( degenes( Tr.Lo.VS.Tr.Na, q= 0.99, M= NULL ) ) ),
                             yes= TRUE, no= FALSE ) )

dim( DE99_biREP.TMM ) #109 genes


# DE triREP (TMM normalized)

DE99_triREP.TMM <- subset( myTMM.triREP,
                           ifelse( row.names( myTMM.triREP ) %in% c(
                             rownames( degenes( Op.Na.VS.Tr.Na, q= 0.99, M= NULL ) ),
                             rownames( degenes( Op.2X.VS.Tr.2X, q= 0.99, M= NULL ) ),
                             rownames( degenes( Op.Na.VS.Op.2X, q= 0.99, M= NULL ) ),
                             rownames( degenes( Tr.Na.VS.Tr.2X, q= 0.99, M= NULL ) ) ),
                             yes= TRUE, no= FALSE ) )

dim( DE99_triREP.TMM ) #190 genes


# ALL (biREP + triREP)

DE99_ALL.TMM <- subset( myTMM,
                        ifelse( row.names( myTMM ) %in% c(
                          rownames( DE99_biREP.TMM ),
                          rownames( DE99_triREP.TMM ) ),
                          yes= TRUE, no= FALSE ) ) 

dim(DE99_ALL.TMM) #249 genes


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


DE99_biREP.TMM.d <- subset( DE99_biREP.TMM,
                            ifelse( row.names( DE99_biREP.TMM ) %in% c(
                              rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel1"]>2.5),] ),
                              rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel2"]>95),] ),
                              rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel3"]>9000),] ),
                              rownames( stats4filt_rslt.biREP.TMM[which( stats4filt_rslt.biREP.TMM[,"rel4"]>200),] )
                            ),
                            yes = FALSE, no = TRUE))

dim(DE99_biREP.TMM.d) #106

#

DE99_triREP.TMM.d <- subset( DE99_triREP.TMM,
                             ifelse( row.names( DE99_triREP.TMM ) %in% c(
                               rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel1"]>2.5),] ),
                               rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel2"]>95),] ),
                               rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel3"]>9000),] ),
                               rownames( stats4filt_rslt.triREP.TMM[which( stats4filt_rslt.triREP.TMM[,"rel4"]>200),] )
                             ),
                             yes = FALSE, no = TRUE))

dim(DE99_triREP.TMM.d) #188

#

DE99_ALL.TMM.d <- subset( myTMM,
                          ifelse( row.names( myTMM ) %in% c(
                            rownames( DE99_biREP.TMM.d ),
                            rownames( DE99_triREP.TMM.d ) ),
                            yes= TRUE, no= FALSE ) ) 

dim(DE99_ALL.TMM.d) #244


# DE biREP - final dataset with VST normalization for plots

DE99_biREP.vst <- subset( assay(vst_biREP.f),
                          ifelse( row.names( assay(vst_biREP.f) ) %in% c(
                            rownames( DE99_biREP.TMM.d ) ),
                            yes= TRUE, no= FALSE ) )

# exploratory PCA for DE biREP 
pcDat.DE99.biREP <- prcomp(t(DE99_biREP.vst))
# plot PCA
autoplot(pcDat.DE99.biREP, data = sample_info.biREP,
         colour="Protein",
         shape="Cactus",
         size=5)

# exploratory heatmap for DE biREP 
hmcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(288)
H<-function(d) hclust(d, method="ward.D2")
heatmap( DE99_biREP.vst[, c(
  "Op.Lo_A","Op.Lo_B","Tr.Lo_A","Tr.Lo_B",
  "Op.Na_A","Op.Na_B","Tr.Na_A","Tr.Na_A" )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )


# DE triREP - final dataset with VST normalization for plots

DE99_triREP.vst <- subset( assay(vst_triREP.f),
                           ifelse( row.names( assay(vst_triREP.f) ) %in% c(
                             rownames( DE99_triREP.TMM.d ) ),
                             yes= TRUE, no= FALSE ) )

# exploratory PCA for DE triREP 
pcDat.DE99.triREP <- prcomp(t(DE99_triREP.vst))
# plot PCA
autoplot(pcDat.DE99.triREP, data = sample_info.triREP,
         colour="Alk.Conc.",
         shape="Cactus",
         size=5)

# exploratory heatmap for DE triREP  
hmcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(288)
H<-function(d) hclust(d, method="ward.D2")
heatmap( DE99_triREP.vst[, c(
  "Op.Na_A","Op.Na_B","Op.Na_C","Tr.Na_A","Tr.Na_B","Tr.Na_C",
  "Op.2X_A","Op.2X_B","Op.2X_C","Tr.2X_A","Tr.2X_B","Tr.2X_C" )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )


# ALL - final dataset with VST normalization for plots

coldata <- data.frame( Sample= Db.sample_info$Sample,
                       Treatment= Db.sample_info$Treatment )

dds <- DESeqDataSetFromMatrix(countData = mycounts_ST,
                              colData = coldata,
                              design = ~ Treatment)
dds

vst.f <- vst(dds, blind=FALSE)
head(assay(vst.f))


DE99_ALL.vst <- subset( assay(vst.f),
                        ifelse( row.names( assay(vst.f) ) %in% c(
                          rownames( DE99_ALL.TMM.d ) ),
                          yes= TRUE, no= FALSE ) )            
dim(DE99_ALL.vst)


hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" ) 
heatmap( DE99_ALL.vst[, c(
  "Op.Lo_A","Op.Lo_B","Tr.Lo_A","Tr.Lo_B",
  "Op.Na_A","Op.Na_B","Op.Na_C","Tr.Na_A","Tr.Na_B","Tr.Na_C",
  "Op.2X_A","Op.2X_B","Op.2X_C","Tr.2X_A","Tr.2X_B","Tr.2X_C" )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

pcDat.DE99.ALL <- prcomp(t(DE99_ALL.vst))
autoplot(pcDat.DE99.ALL, data = sample_info,
         colour="Alk.Conc.",
         shape="Cactus",
         size=5) 



# MDS ##################################################

d_ALL <- dist( t( DE99_ALL.vst ) )
MDS_ALL <- smacofSym(d_ALL, ndim=3) #metric iterative MDS, no PCA
MDS_ALL$stress #Low stress, the analysis is OK

rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = "white" )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.8)
}

rgl_init()

rgl.spheres( MDS_ALL$conf[1:3,1], MDS_ALL$conf[1:3,2], MDS_ALL$conf[1:3,3],
             color= "#9e0142", r= 0.025, lit= F ) #Op.2X_Dbuz
rgl.spheres( MDS_ALL$conf[4:6,1], MDS_ALL$conf[4:6,2], MDS_ALL$conf[4:6,3],
             color= "#f46d43", r= 0.025, lit= F ) #Op.Na_Dbuz
rgl.spheres( MDS_ALL$conf[7:8,1], MDS_ALL$conf[7:8,2], MDS_ALL$conf[7:8,3],
             color= "#ffc762", r= 0.025, lit= F ) #Op.Lo_Dbuz

rgl.spheres( MDS_ALL$conf[9:11,1], MDS_ALL$conf[9:11,2], MDS_ALL$conf[9:11,3],
             color= "#3d85c6", r= 0.025, lit= F ) #Tr.2X_Dbuz
rgl.spheres( MDS_ALL$conf[12:14,1], MDS_ALL$conf[12:14,2], MDS_ALL$conf[12:14,3],
             color= "#66c2a5", r= 0.025, lit= F ) #Tr.Na_Dbuz
rgl.spheres( MDS_ALL$conf[15:16,1], MDS_ALL$conf[15:16,2], MDS_ALL$conf[15:16,3],
             color= "#e6f598", r= 0.025, lit= F ) #Tr.Lo_Dbuz


triangles3d( MDS_ALL$conf[1:3,1], MDS_ALL$conf[1:3,2], MDS_ALL$conf[1:3,3],
            color= "#9e0142", lit=F, alpha=0.25, lwd=2)
triangles3d( MDS_ALL$conf[4:6,1], MDS_ALL$conf[4:6,2], MDS_ALL$conf[4:6,3],
            color= "#f46d43", lit=F, alpha=0.25, lwd=2)
segments3d( MDS_ALL$conf[7:8,1], MDS_ALL$conf[7:8,2], MDS_ALL$conf[7:8,3],
            color= "#ffc762", lit=F, alpha=0.3, lwd=3)

triangles3d( MDS_ALL$conf[9:11,1], MDS_ALL$conf[9:11,2], MDS_ALL$conf[9:11,3],
            color= "#3d85c6", lit=F, alpha=0.25, lwd=2)
triangles3d( MDS_ALL$conf[12:14,1], MDS_ALL$conf[12:14,2], MDS_ALL$conf[12:14,3],
            color= "#66c2a5", lit=F, alpha=0.25, lwd=2)
segments3d( MDS_ALL$conf[15:16,1], MDS_ALL$conf[15:16,2], MDS_ALL$conf[15:16,3],
            color= "#e6f598", lit=F, alpha=0.3, lwd=3)


ellipse_0X <- ellipse3d( cov( cbind( MDS_ALL$conf[c( 4:6,7:8 ), 1],
                                       MDS_ALL$conf[c( 4:6,7:8 ), 2],
                                       MDS_ALL$conf[c( 4:6,7:8 ), 3]) ),
                           centre= c(
                             mean( MDS_ALL$conf[c( 4:6,7:8 ), 1]),
                             mean( MDS_ALL$conf[c( 4:6,7:8 ), 2]),
                             mean( MDS_ALL$conf[c( 4:6,7:8 ), 3]) ),
                           scale= c( 1,1,1 ), level= 0.95 )

wire3d( ellipse_0X, col= "#FA9A53",  lit= F, alpha= 0.1 ) #0X alkaloids 

ellipse_1X <- ellipse3d( cov( cbind( MDS_ALL$conf[c( 12:14,15:16 ), 1],
                                     MDS_ALL$conf[c( 12:14,15:16 ), 2],
                                     MDS_ALL$conf[c( 12:14,15:16 ), 3]) ),
                         centre= c(
                           mean( MDS_ALL$conf[c( 12:14,15:16 ), 1]),
                           mean( MDS_ALL$conf[c( 12:14,15:16 ), 2]),
                           mean( MDS_ALL$conf[c( 12:14,15:16 ), 3]) ),
                         scale= c( 1,1,1 ), level= 0.95 )

wire3d( ellipse_1X, col= "#A6DC9F",  lit= F, alpha= 0.1 ) #1X alkaloids

ellipse_2X <- ellipse3d( cov( cbind( MDS_ALL$conf[c( 1:3,9:11 ), 1],
                                     MDS_ALL$conf[c( 1:3,9:11 ), 2],
                                     MDS_ALL$conf[c( 1:3,9:11 ), 3]) ),
                         centre= c(
                           mean( MDS_ALL$conf[c( 1:3,9:11 ), 1]),
                           mean( MDS_ALL$conf[c( 1:3,9:11 ), 2]),
                           mean( MDS_ALL$conf[c( 1:3,9:11 ), 3]) ),
                         scale= c( 1,1,1 ), level= 0.95 )

wire3d( ellipse_2X, col= "#6E4384",  lit= F, alpha= 0.1 ) #2X alkaloids 

axes3d( edges= "bbox", labels= F, tick= F, box= F, expand= 1.025, color= "#252525" )
rgl.material(smooth = T, point_antialias = T, line_antialias = T)


#TOP100 most variable DE genes

hmcol <- colorRampPalette( brewer.pal( 9, "BuPu" ) )( 288 )
H <- function( d ) hclust( d, method= "ward.D2" ) 

# Estimation of variance for each row in the log-counts matrix
countVar_ALL <- apply(DE99_ALL.vst, 1, var)
# Get row numbers for the top 100 most variable
highVar_ALL <- order(countVar_ALL, decreasing=TRUE)[1:100]
# Subset log-counts matrix
top100_ALL <- DE99_ALL.vst[highVar_ALL,]

hm.TOP100 <- heatmap( top100_ALL[, c(
  "Op.Lo_A","Op.Lo_B","Tr.Lo_A","Tr.Lo_B",
  "Op.Na_A","Op.Na_B","Op.Na_C","Tr.Na_A","Tr.Na_B","Tr.Na_C",
  "Op.2X_A","Op.2X_B","Op.2X_C","Tr.2X_A","Tr.2X_B","Tr.2X_C"  )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( top100_ALL )[ hm.TOP100$rowInd ] ),
             "Db.hm.TOP100_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


# Heatmap of DE genes from enzime categories
# Import DE genes from evaluated enzime categories
ez <- read.table(file = "Db_DEenzyme.names", header = F, row.names = 1)

ez_vst <- subset( DE99_ALL.vst,
                    ifelse( row.names( DE99_ALL.vst ) %in% c( rownames( ez ) ),
                            yes= TRUE, no= FALSE ) )

hm.ez <- heatmap( ez_vst[, c(
  "Op.Lo_A","Op.Lo_B","Tr.Lo_A","Tr.Lo_B",
  "Op.Na_A","Op.Na_B","Op.Na_C","Tr.Na_A","Tr.Na_B","Tr.Na_C",
  "Op.2X_A","Op.2X_B","Op.2X_C","Tr.2X_A","Tr.2X_B","Tr.2X_C"  )],
  Colv= NA, col= hmcol, cexRow= 0.5, cexCol= 1, hclustfun= H )

write.table( rev( rownames( ez_vst )[ hm.ez$rowInd ] ),
             "Db.hm.ez_rows",
             quote= FALSE, sep= "\t", row.names= FALSE, col.names= FALSE )


