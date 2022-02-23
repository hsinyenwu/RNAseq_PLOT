
### `PLOTmp` is a function to plot multiple RNAseq datasets.

It requires the following bioconductor libraries:
```
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
```

Example of using `gene.structure` for loading a gff file (also can use gtf files, but need to change `format="gtf`")

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
(Here I only have three files, but use them five times to show the function can take more files.
The number of RNA-seq files you can plot simultaneously depend on the RAM size of your computer.)
```
RNAseqData1 = cbind(CTRL_selected,ABA_all,CTRL_all,ABA_all,CTRL_selected)
```

```
PLOTmp(AT1G01740,RNAseqData=RNAseqData1,Extend=200,BGcolor="#FEFEAE")
```
<img width="739" alt="image" src="https://user-images.githubusercontent.com/4383665/153487394-df0cfd86-46fb-4fe8-8ca8-c0a77e06aee0.png">

Parameters for `PLOTmp`
```
#' @YFG: gene id 
#' @RNAseqData: the vector of the bam file paths
#' @Extend x number of nucleotide on both side of YFG (wider on both sides)
#' @BGcolor: RNA-seq coverage color

```

```
sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
 [1] rtracklayer_1.50.0          GenomicAlignments_1.26.0   
 [3] Rsamtools_2.6.0             Biostrings_2.58.0          
 [5] XVector_0.30.0              SummarizedExperiment_1.20.0
 [7] MatrixGenerics_1.2.1        matrixStats_0.60.1         
 [9] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0       
[11] Biobase_2.50.0              GenomicRanges_1.42.0       
[13] GenomeInfoDb_1.26.7         IRanges_2.24.1             
[15] S4Vectors_0.28.1            BiocGenerics_0.36.1        

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7             lattice_0.20-44        prettyunits_1.1.1     
 [4] assertthat_0.2.1       utf8_1.2.2             BiocFileCache_1.14.0  
 [7] R6_2.5.1               RSQLite_2.2.8          httr_1.4.2            
[10] pillar_1.6.2           zlibbioc_1.36.0        rlang_0.4.11          
[13] progress_1.2.2         curl_4.3.2             rstudioapi_0.13       
[16] blob_1.2.2             Matrix_1.3-4           BiocParallel_1.24.1   
[19] stringr_1.4.0          RCurl_1.98-1.4         bit_4.0.4             
[22] biomaRt_2.46.3         DelayedArray_0.16.3    compiler_4.0.3        
[25] pkgconfig_2.0.3        askpass_1.1            openssl_1.4.5         
[28] tidyselect_1.1.1       tibble_3.1.4           GenomeInfoDbData_1.2.4
[31] XML_3.99-0.7           fansi_0.5.0            crayon_1.4.1          
[34] dplyr_1.0.7            dbplyr_2.1.1           bitops_1.0-7          
[37] rappdirs_0.3.3         grid_4.0.3             lifecycle_1.0.0       
[40] DBI_1.1.1              magrittr_2.0.1         stringi_1.7.4         
[43] cachem_1.0.6           xml2_1.3.2             ellipsis_0.3.2        
[46] generics_0.1.0         vctrs_0.3.8            tools_4.0.3           
[49] bit64_4.0.5            glue_1.4.2             purrr_0.3.4           
[52] hms_1.1.0              yaml_2.2.1             fastmap_1.1.0         
[55] memoise_2.0.0      
```
