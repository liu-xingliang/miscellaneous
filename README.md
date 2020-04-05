# install monocle3 and seurat3 (develop) under conda

Create virtual environment in anaconda

```bash
conda create -n R363 anaconda
```

Make sure `$LD_LIBRARY_PATH` is correctly setted when virual env is activated and deactivated 

```bash
cd /home/liuxl18/anaconda3/envs/R363
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
```

`activate.d/env_vars.sh`:
```bash
export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/home/liuxl18/anaconda3/envs/R363/lib:/home/liuxl18/anaconda3/lib:${LD_LIBRARY_PATH}
```

`deactivate.d/env_vars.sh`:
```bash
export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}
unset OLD_LD_LIBRARY_PATH
```

Install necessary packages for Monocle3 and Seurat3 (hdf5)

```bash
conda activate R363
conda install -c conda-forge r-base=3.6.3 gdal proj udunits
conda install -c anaconda geos hdf5
ln -s /home/liuxl18/anaconda3/envs/R363/lib/libzstd.so.1.3.7 /home/liuxl18/anaconda3/envs/R363/lib/libzstd.so.1
```

Running in R

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
                       
# devtools
install.packages("devtools", repos="https://cloud.r-project.org")

# monocle3
library(devtools)
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

# seurat3 (develop)
devtools::install_github("satijalab/seurat", ref = "develop")
```

`sessionInfo()`:

```r
library(monocle3)
library(Seurat)

> library(monocle3)
Loading required package: Biobase
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:parallel’:

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unsplit, which, which.max, which.min

Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.

Loading required package: SingleCellExperiment
Loading required package: SummarizedExperiment
Loading required package: GenomicRanges
Loading required package: stats4

Loading required package: S4Vectors

Attaching package: ‘S4Vectors’

The following object is masked from ‘package:base’:

    expand.grid

Loading required package: IRanges
Loading required package: GenomeInfoDb
Loading required package: DelayedArray
Loading required package: matrixStats

Attaching package: ‘matrixStats’

The following objects are masked from ‘package:Biobase’:

    anyMissing, rowMedians

Loading required package: BiocParallel

Attaching package: ‘DelayedArray’

The following objects are masked from ‘package:matrixStats’:

    colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

The following objects are masked from ‘package:base’:

    aperm, apply, rowsum


Attaching package: ‘monocle3’

The following objects are masked from ‘package:Biobase’:

    exprs, fData, fData<-, pData, pData<-

> 
> library(Seurat)

Attaching package: ‘Seurat’

The following object is masked from ‘package:SummarizedExperiment’:

    Assays

sessionInfo()

> sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /home/liuxl18/anaconda3/envs/R363/lib/libmkl_rt.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] Seurat_3.1.4.9020           monocle3_0.2.1             
 [3] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
 [5] DelayedArray_0.12.2         BiocParallel_1.20.1        
 [7] matrixStats_0.56.0          GenomicRanges_1.38.0       
 [9] GenomeInfoDb_1.22.1         IRanges_2.20.2             
[11] S4Vectors_0.24.3            Biobase_2.46.0             
[13] BiocGenerics_0.32.0         devtools_2.2.2             
[15] usethis_1.5.1              

loaded via a namespace (and not attached):
  [1] Rtsne_0.15             colorspace_1.4-1       ellipsis_0.3.0        
  [4] ggridges_0.5.2         rprojroot_1.3-2        XVector_0.26.0        
  [7] fs_1.4.1               leiden_0.3.3           listenv_0.8.0         
 [10] npsurv_0.4-0           remotes_2.1.1          ggrepel_0.8.2         
 [13] fansi_0.4.1            codetools_0.2-16       splines_3.6.3         
 [16] lsei_1.2-0             pkgload_1.0.2          jsonlite_1.6.1        
 [19] ica_1.0-2              cluster_2.1.0          png_0.1-7             
 [22] uwot_0.1.8             sctransform_0.2.1      BiocManager_1.30.10   
 [25] compiler_3.6.3         httr_1.4.1             backports_1.1.6       
 [28] lazyeval_0.2.2         assertthat_0.2.1       Matrix_1.2-18         
 [31] cli_2.0.2              htmltools_0.4.0        prettyunits_1.1.1     
 [34] tools_3.6.3            rsvd_1.0.3             igraph_1.2.5          
 [37] gtable_0.3.0           glue_1.4.0             GenomeInfoDbData_1.2.2
 [40] RANN_2.6.1             reshape2_1.4.3         dplyr_0.8.5           
 [43] rappdirs_0.3.1         Rcpp_1.0.4             vctrs_0.2.4           
 [46] gdata_2.18.0           ape_5.3                nlme_3.1-145          
 [49] lmtest_0.9-37          stringr_1.4.0          globals_0.12.5        
 [52] ps_1.3.2               testthat_2.3.2         lifecycle_0.2.0       
 [55] irlba_2.3.3            gtools_3.8.2           future_1.16.0         
 [58] zlibbioc_1.32.0        MASS_7.3-51.5          zoo_1.8-7             
 [61] scales_1.1.0           RColorBrewer_1.1-2     curl_4.3              
 [64] memoise_1.1.0          reticulate_1.15        pbapply_1.4-2         
 [67] gridExtra_2.3          ggplot2_3.3.0          stringi_1.4.6         
 [70] desc_1.2.0             caTools_1.18.0         pkgbuild_1.0.6        
 [73] rlang_0.4.5            pkgconfig_2.0.3        bitops_1.0-6          
 [76] lattice_0.20-41        ROCR_1.0-7             purrr_0.3.3           
 [79] htmlwidgets_1.5.1      patchwork_1.0.0        cowplot_1.0.0         
 [82] processx_3.4.2         tidyselect_1.0.0       RcppAnnoy_0.0.16      
 [85] plyr_1.8.6             magrittr_1.5           R6_2.4.1              
 [88] gplots_3.0.3           pillar_1.4.3           withr_2.1.2           
 [91] fitdistrplus_1.0-14    survival_3.1-11        RCurl_1.98-1.1        
 [94] tsne_0.1-3             tibble_3.0.0           future.apply_1.4.0    
 [97] crayon_1.3.4           KernSmooth_2.23-16     plotly_4.9.2.1        
[100] viridis_0.5.1          grid_3.6.3             data.table_1.12.8     
[103] callr_3.4.3            digest_0.6.25          tidyr_1.0.2           
[106] munsell_0.5.0          viridisLite_0.3.0      sessioninfo_1.1.1    
```


