```
###############################################
#rm(list=ls())

##################################################
# Only plot the riboseq reads in the CDS region  #
##################################################

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)

# gene.structure function help to load all required information to the global environment
# plotRanges plots GRanges for an isoform 
# plotGeneModel_num uses plotRanges to plot all isoforms
# plot multiple paired-end data

##################################################
##################################################
##################################################
##################################################
gene.structure <- function(annotation,format="gtf",dataSource="",organism=""){
  txdb <- makeTxDbFromGFF(annotation,format=format, dataSource=dataSource,organism=organism)
  exonsByTx <- exonsBy(txdb,by='tx',use.names=T)
  exonsByGene <- exonsBy(txdb,by='gene')
  exonsByTxGene <- exonsBy(txdb,by=c('tx','gene'),use.names=T)
  txByGene <- transcriptsBy(txdb,by='gene')
  cdsByTx <- cdsBy(txdb, by="tx",use.names=T)
  fiveUTR <- fiveUTRsByTranscript(txdb,use.names=T)
  threeUTR <- threeUTRsByTranscript(txdb,use.names=T)
  cds <- cdsBy(txdb, by=c("tx","gene"),use.names=TRUE)
  
  assign("txdb", txdb, envir = .GlobalEnv)
  assign("exonsByTx", exonsByTx, envir = .GlobalEnv)
  assign("exonsByGene", exonsByGene, envir = .GlobalEnv)
  assign("exonsByTxGene", exonsByTxGene, envir = .GlobalEnv)
  assign("txByGene", txByGene, envir = .GlobalEnv)
  assign("cdsByTx", cdsByTx, envir = .GlobalEnv)
  assign("fiveUTR", fiveUTR, envir = .GlobalEnv)
  assign("threeUTR", threeUTR, envir = .GlobalEnv)
  assign("cds",cds, envir = .GlobalEnv)
}

plotRanges <- function(isoform,shortest3UTR, ybottom, main = deparse(substitute(x)),colCDS = "black",col3="white",col5="lightgrey") {
  if(isoform %in% names(cdsByTx)) {
    height <- 0.1
    xlim=ranges(unlist(exonsByTx[isoform]))
    xlimCds=ranges(unlist(cdsByTx[isoform]))
    # plot 5'UTR
    if (isoform %in% names(fiveUTR)) {
      xlim5=ranges(unlist(fiveUTR[isoform]))
      rect(start(xlim5), ybottom, end(xlim5), ybottom + height, col = col5, border = "black")
    }
    # plot lines between exons to represent introns
    if (length(unlist(exonsByTx[isoform]))>1) {
      GAPS <- gaps(unlist(exonsByTx[isoform]),start=NA)
      segments(x0 = start(GAPS),
               y0 = ybottom+height/2,
               x1 = start(GAPS)+width(ranges(GAPS))/2,
               y1 = ybottom+height,
               col = "black",lwd=1)
      
      segments(x0 = start(GAPS)+width(ranges(GAPS))/2,
               y0 = ybottom+height,
               x1 = end(GAPS),
               y1 = ybottom+height/2,
               col = "black",lwd=1)
    }
    
    rect(start(xlimCds), ybottom, end(xlimCds), ybottom + height, col =c("black","black","black") , border = "black")
    
    # Plot 3'UTR with an arrow shap
    if (isoform %in% names(threeUTR)) {
      xlim3=ranges(sort(unlist(threeUTR[isoform])))
      Length=length(unlist(threeUTR[isoform]))
      if (shortest3UTR <=50) {
        z=shortest3UTR
      }
      else{
        z=shortest3UTR/3
      }
      if (as.character(runValue(strand(exonsByTx[isoform])))=="+") {
        if (length(unlist(threeUTR[isoform]))==1) {
          polygon(x=c(start(xlim3), end(xlim3)-z,end(xlim3),end(xlim3)-z,start(xlim3)),y=c(ybottom+height,ybottom+height,ybottom+height/2,ybottom,ybottom), col = col3, border = "black")
        }
        else {
          rect(start(xlim3[1:Length-1]), ybottom, end(xlim3[1:Length-1]), ybottom + height, col = col3, border = "black")
          polygon(x=c(start(xlim3[Length]), end(xlim3[Length])-z,end(xlim3[Length]),end(xlim3[Length])-z,start(xlim3[Length])),y=c(ybottom+height,ybottom+height,ybottom+height/2,ybottom,ybottom), col = col3, border = "black")
        }
      }
      if (as.character(runValue(strand(exonsByTx[isoform])))=="-") {
        if (length(unlist(threeUTR[isoform]))==1) {
          polygon(x=c(start(xlim3),start(xlim3)+z,end(xlim3),end(xlim3),start(xlim3)+z),y=c(ybottom+height/2,ybottom+height,ybottom+height,ybottom,ybottom), col = col3, border = "black")
        }
        else {
          rect(start(xlim3[2:Length]), ybottom, end(xlim3[2:Length]), ybottom+height, col = col3, border = "black")
          polygon(x=c(start(xlim3[1]),start(xlim3[1])+z,end(xlim3[1]),end(xlim3[1]),start(xlim3[1])+z),y=c(ybottom+height/2,ybottom+height,ybottom+height,ybottom,ybottom), col = col3, border = "black")
        }
      }
    }
    axis(1)
  }
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
}

plotGeneModel_num <- function(gene,Extend=Extend){
  isoforms <- length(unlist(txByGene[gene]))
  generanges <- ranges(unlist(exonsByGene[gene]))
  SUW <- sum(width(generanges))
  xlimg= min(start(generanges))-0.05
  genelim <- c(min(start(generanges))-Extend, max(end(generanges))+Extend)
  isoforms.w.3UTR <- unlist(txByGene[gene])$tx_name[which(unlist(txByGene[gene])$tx_name %in% names(threeUTR))]
  plot.new()
  yAxis <- (isoforms*0.3+0.1)
  plot.window(genelim,c(0,yAxis))
  tx_name_start_pos <- nchar(gene)+2 #find the position of the tx name start
  tx_num <- sort(substr(unlist(txByGene[gene])$tx_name,tx_name_start_pos,nchar(unlist(txByGene[gene])$tx_name)))
  # tx_fac <- as.numeric(as.factor(tx_num))
  for (i in sort(unlist(txByGene[gene])$tx_name)) {
    k=as.numeric(substr(i,tx_name_start_pos,nchar(i)))
    k2=which(tx_num==k)
    if (i %in% names(threeUTR)) {
      shortest3UTR <- min(sapply(isoforms.w.3UTR, function(j) width(tail(unlist(threeUTR[j]),1))))
      plotRanges(isoform=i,shortest3UTR,ybottom=(yAxis-0.28*k2)) #removed
      text(x=min(start(generanges))-Extend-0.1, y=(yAxis-0.28*k2+0.05), labels=tx_num[k2],cex=1.2)
    }
    else {
      plotRanges(isoform=i,ybottom=(yAxis-0.28*k2))
      text(x=min(start(generanges))-Extend-0.1, y=(yAxis-0.28*k2+0.05), labels=tx_num[k2],cex=1.2)
    }
  }
}

PLOTmp <-function(YFG,RNAseqData,Extend=0,BGcolor="#FEFEAE") {
  #RNAseqData is a list of defined path (e.g. RNAseqData-list(RNAseqBam1,RNAseqBam2,RNAseqBam3)
  #You have to define RNAseqBam1 <- "~/Desktop/RNA1.bam" ...
  YFG <- deparse(substitute(YFG)) #So you do not need to add quote "" around the gene name
  Names_Dataset <- colnames(RNAseqData)
  Num_Dataset <- length(RNAseqData)
  Num_row <- Num_Dataset+1
  
  par(mfrow=c(Num_row,1),mar=c(0.2,0.2,0.2,0.2),oma=c(3,2,4,2))
  
  chr <- as.character(seqnames(exonsByGene[YFG])[[1]])[1]
  generanges <- ranges(unlist(exonsByGene[YFG]))
  GR <- GRanges(seqnames=as.character(chr),IRanges(min(start(generanges))-Extend, max(end(generanges))+Extend),strand=strand(unlist(exonsByGene[YFG]))[1])
  #Add ranges for extracting RNAseq reads
  SZ <- GenomicRanges::reduce(unlist(txByGene[YFG]))
  which1 <- resize(SZ,width=width(SZ)+Extend,fix = "end")
  which1 <- resize(which1,width=width(which1)+Extend,fix = "start")
  what1 <- c("rname", "strand", "pos", "qwidth","seq")
  param <- ScanBamParam(which = which1, what = what1)
  
  isoforms <- length(unlist(txByGene[YFG]))
  layout(matrix(rep(seq(1:Num_row),each=2),Num_row,2,byrow=TRUE), widths=c(rep(6,Num_row)), heights=c(rep(3,Num_Dataset),(0.4+Num_row*0.08)*isoforms))
  
  for(i in seq_len(Num_Dataset)){
    readPairs <- readGAlignmentPairs(RNAseqData[i], param=param,strandMode = 2)
    readPairs <- readPairs[strand(readPairs)==as.character(strand(GR))]
    cvg <- coverage(readPairs)
    Gtx <- as.numeric(cvg[[chr]][ranges(GR)])
    plot(Gtx,type="h",col=BGcolor,lwd=1,xaxt='n',ylim=c(0,max(Gtx)+2))
    # mtext(RNAlab1, side = 2, line = 2.5,cex=1)
    par(new = T)
    plot(Gtx,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,max(Gtx)+2))
    lines(x=c(1,length(Gtx)),y=c(0,0),col="white",lwd=2)
    legend("topleft",inset=c(-0.03,-0.03),legend=Names_Dataset[i],bty="n",cex=1.5)
  }
  plotGeneModel_num(YFG,Extend=Extend)
  mtext(YFG,side=3,line=0.4, cex=1.2, col="black", outer=TRUE)
}

##################################################
##################################################
##################################################
##################################################

# The following show you how to run the code
############################
#  Load RNASeq Bam files   #
############################

# Variables for RNA-seq bam file paths 
# The names of the variables will be shown on the topleft of each RNA-seq coverage plot 

CTRL_selected="~/Desktop/ABA/P_sites_all/RNA_DMSO60_merged_selected_0.935.bam"
ABA_all="~/Desktop/ABA/P_sites_all/RNA_ABA60_merged.bam"
CTRL_all="~/Desktop/CTRL_v1/RNA_CTRL_merged.bam"

# make a vector for the paths of the bam files, this 
RNAseqData1 = cbind(CTRL_selected,ABA_all,CTRL_all,ABA_all,CTRL_selected)

############################
#  Start plotting genes   #
############################

# gene.structure(annotation="~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
gene.structure(annotation="~/Desktop/Arabidopsis/Araport11_GFF3_genes_transposons.201606.gff",
               format="gff",
               dataSource="Araport",
               organism="Arabidopsis thaliana")
# PLOTmp parameters
#@ YFG: gene id 
#@ RNAseqData: the vector of the bam file paths
#@ Extend x number of nucleotide on both side of YFG (wider on both sides)
#@ BGcolor: RNA-seq coverage color

PLOTmp(YFG=AT1G01060,RNAseqData=RNAseqData1,Extend=100,BGcolor="#FEFEAE")
PLOTmp(AT1G01020,RNAseqData=RNAseqData1,Extend=2000,BGcolor="lightgrey")
PLOTmp(AT1G01370,RNAseqData=RNAseqData1,Extend=0)
PLOTmp(AT1G01740,RNAseqData=RNAseqData1,Extend=200)
PLOTmp(AT1G01220,RNAseqData=RNAseqData1,Extend=200)
PLOTmp(AT1G01700,RNAseqData=RNAseqData1,Extend=200)
```
