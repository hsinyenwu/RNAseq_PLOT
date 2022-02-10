### RNAseq_PLOT
PLOTmp is a function to plot any number of RNAseq datasets.

It requires the following libraries from bioconductor:
```
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
```

```
#Example gff loading (also can use gtf files, but need to change format="gtf"

gene.structure(annotation="~/Desktop/Arabidopsis/Araport11_GFF3_genes_transposons.201606.gff",
               format="gff",
               dataSource="Araport",
               organism="Arabidopsis thaliana")
```

```
PLOTmp(AT1G01740,RNAseqData=RNAseqData1,Extend=200,BGcolor="#FEFEAE")
```
<img width="739" alt="image" src="https://user-images.githubusercontent.com/4383665/153487394-df0cfd86-46fb-4fe8-8ca8-c0a77e06aee0.png">
