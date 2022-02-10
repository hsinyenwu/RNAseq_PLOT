### RNAseq_PLOT
PLOTmp is a function to plot any number of RNAseq datasets.

It requires the following libraries from bioconductor:
```
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
```

Example of using `gene.structure` for loading a gff file (also can use gtf files, but need to change format="gtf"

```
gene.structure(annotation="~/Desktop/Arabidopsis/Araport11_GFF3_genes_transposons.201606.gff",
               format="gff",
               dataSource="Araport",
               organism="Arabidopsis thaliana")
```

Next, create variables for RNA-seq bam file paths 
```
# The names of the variables will be shown on the topleft of each RNA-seq coverage plot 

CTRL_selected="~/Desktop/ABA/P_sites_all/RNA_DMSO60_merged_selected_0.935.bam"
ABA_all="~/Desktop/ABA/P_sites_all/RNA_ABA60_merged.bam"
CTRL_all="~/Desktop/CTRL_v1/RNA_CTRL_merged.bam"
```

Next, make a vector for the paths of the bam files
```
# Here I only have three files, but use them five times to show the function can take more files.
# How many RNA-seq files you can plot simultaneously depend on the RAM size of your computer.
RNAseqData1 = cbind(CTRL_selected,ABA_all,CTRL_all,ABA_all,CTRL_selected)
```

```
PLOTmp(AT1G01740,RNAseqData=RNAseqData1,Extend=200,BGcolor="#FEFEAE")
```
<img width="739" alt="image" src="https://user-images.githubusercontent.com/4383665/153487394-df0cfd86-46fb-4fe8-8ca8-c0a77e06aee0.png">
