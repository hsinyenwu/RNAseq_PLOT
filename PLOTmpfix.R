#Same as PLOTmp except, for one given gene, use the max of the first RNAseq dataset for all datasets

PLOTmpfix <-function(YFG,RNAseqData,Extend=0,BGcolor="#FEFEAE") {
  #RNAseqData is a list of defined path (e.g. RNAseqData-list(RNAseqBam1,RNAseqBam2,RNAseqBam3)
  #You have to define RNAseqBam1 <- "~/Desktop/RNA1.bam" ...
  # YFG <- deparse(substitute(YFG)) #So you do not need to add quote "" around the gene name
  print(paste("gene name:",YFG))
  Names_Dataset <- colnames(RNAseqData)
  Num_Dataset <- length(RNAseqData)
  Num_row <- Num_Dataset+1
  
  par(mfrow=c(Num_row,1),mar=c(0.2,0.2,0.2,0.2),oma=c(3,2,4,2))
  
  chr <- as.character(seqnames(exonsByGene[YFG])[[1]])[1]
  # print(chr)
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
    if(i==1){Gtx_n <- max(Gtx)}
    plot(Gtx,type="h",col=BGcolor,lwd=1,xaxt='n',ylim=c(0,Gtx_n+2))
    # mtext(RNAlab1, side = 2, line = 2.5,cex=1)
    par(new = T)
    plot(Gtx,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,Gtx_n+2))
    lines(x=c(1,length(Gtx)),y=c(0,0),col="white",lwd=2)
    legend("topleft",legend=Names_Dataset[i],bty="n",cex=1) #inset=c(-0.05,-0.05)
  }
  plotGeneModel_num(YFG,Extend=Extend)
  mtext(YFG,side=3,line=0.4, cex=1.2, col="black", outer=TRUE)
}
