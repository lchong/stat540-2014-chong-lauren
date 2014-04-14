Seminar 7 -- RNA-seq: differential expression analysis
======================================================
### Lauren Chong

Looking at read count data from RNA-seq of 21 mice (2 strains).

Load packages:

```r
library(edgeR)
library(DESeq)
library(limma)
```


Load the data:

```r
dat <- read.table("bottomly_count_table.tsv", header = TRUE, row.names = 1)
des <- read.table("bottomly_phenodata.tsv", header = TRUE, row.names = 1)
dim(dat)
```

```
## [1] 36536    21
```

```r
head(dat)
```

```
##                    SRX033480 SRX033488 SRX033481 SRX033489 SRX033482
## ENSMUSG00000000001       369       744       287       769       348
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         0         1         0         1         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         0         1         1         5         0
## ENSMUSG00000000049         0         1         0         1         0
##                    SRX033490 SRX033483 SRX033476 SRX033478 SRX033479
## ENSMUSG00000000001       803       433       469       585       321
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         1         0         7         6         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         4         0         0         0         0
## ENSMUSG00000000049         0         0         0         0         0
##                    SRX033472 SRX033473 SRX033474 SRX033475 SRX033491
## ENSMUSG00000000001       301       461       309       374       781
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         1         1         1         1         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         4         1         1         0         1
## ENSMUSG00000000049         0         0         0         0         0
##                    SRX033484 SRX033492 SRX033485 SRX033493 SRX033486
## ENSMUSG00000000001       555       820       294       758       419
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         2         1         1         4         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         2         1         1         1         1
## ENSMUSG00000000049         0         0         0         0         0
##                    SRX033494
## ENSMUSG00000000001       857
## ENSMUSG00000000003         0
## ENSMUSG00000000028         5
## ENSMUSG00000000031         0
## ENSMUSG00000000037         2
## ENSMUSG00000000049         0
```

```r
show(des)
```

```
##           num.tech.reps   strain experiment.number lane.number
## SRX033480             1 C57BL/6J                 6           1
## SRX033488             1 C57BL/6J                 7           1
## SRX033481             1 C57BL/6J                 6           2
## SRX033489             1 C57BL/6J                 7           2
## SRX033482             1 C57BL/6J                 6           3
## SRX033490             1 C57BL/6J                 7           3
## SRX033483             1 C57BL/6J                 6           5
## SRX033476             1 C57BL/6J                 4           6
## SRX033478             1 C57BL/6J                 4           7
## SRX033479             1 C57BL/6J                 4           8
## SRX033472             1   DBA/2J                 4           1
## SRX033473             1   DBA/2J                 4           2
## SRX033474             1   DBA/2J                 4           3
## SRX033475             1   DBA/2J                 4           5
## SRX033491             1   DBA/2J                 7           5
## SRX033484             1   DBA/2J                 6           6
## SRX033492             1   DBA/2J                 7           6
## SRX033485             1   DBA/2J                 6           7
## SRX033493             1   DBA/2J                 7           7
## SRX033486             1   DBA/2J                 6           8
## SRX033494             1   DBA/2J                 7           8
```


edgeR analysis
--------------
First we will run our DE analysis using edgeR to model the expression.  

Create a DGEList separating the samples by strain, and then a design matrix for model fitting:

```r
group <- factor(as.numeric(des$strain))
dge.glm <- DGEList(counts = dat, group = group)
show(dge.glm$samples)
```

```
##           group lib.size norm.factors
## SRX033480     1  3040296            1
## SRX033488     1  6303665            1
## SRX033481     1  2717092            1
## SRX033489     1  6545795            1
## SRX033482     1  3016179            1
## SRX033490     1  7097379            1
## SRX033483     1  3707895            1
## SRX033476     1  5165144            1
## SRX033478     1  4953201            1
## SRX033479     1  4192872            1
## SRX033472     2  3970729            1
## SRX033473     2  4733003            1
## SRX033474     2  3702051            1
## SRX033475     2  3569483            1
## SRX033491     2  7276198            1
## SRX033484     2  4422272            1
## SRX033492     2  7115851            1
## SRX033485     2  3467730            1
## SRX033493     2  7339817            1
## SRX033486     2  3879114            1
## SRX033494     2  6771680            1
```

```r
design <- model.matrix(~group)
show(design)
```

```
##    (Intercept) group2
## 1            1      0
## 2            1      0
## 3            1      0
## 4            1      0
## 5            1      0
## 6            1      0
## 7            1      0
## 8            1      0
## 9            1      0
## 10           1      0
## 11           1      1
## 12           1      1
## 13           1      1
## 14           1      1
## 15           1      1
## 16           1      1
## 17           1      1
## 18           1      1
## 19           1      1
## 20           1      1
## 21           1      1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$group
## [1] "contr.treatment"
```


Estimate dispersion for each probe:

```r
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
```

```
## Disp = 0.03893 , BCV = 0.1973
```

```r
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp, design)
```

```
## Loading required package: splines
```

```r
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)
plotBCV(dge.glm.tag.disp)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


Fit an adjusted linear model, run likelihood test:

```r
fit <- glmFit(dge.glm.tag.disp, design)
lrt <- glmLRT(fit, coef = 2)
```


We want to identify genes that show significantly different expression between the two strains. Look at the top 10 most significant results:

```r
topTags(lrt)
```

```
## Coefficient:  group2 
##                     logFC logCPM    LR    PValue       FDR
## ENSMUSG00000020912 -5.187  3.155 401.8 2.238e-89 8.177e-85
## ENSMUSG00000050141 -5.363  2.319 311.6 9.787e-70 1.788e-65
## ENSMUSG00000035775 -4.543  2.674 298.0 9.042e-67 1.101e-62
## ENSMUSG00000015484 -1.968  4.307 282.5 2.113e-63 1.930e-59
## ENSMUSG00000024248 -3.152  3.463 277.0 3.337e-62 2.438e-58
## ENSMUSG00000030532  1.547  5.572 269.4 1.568e-60 9.548e-57
## ENSMUSG00000054354 -6.283  1.940 249.3 3.657e-56 1.909e-52
## ENSMUSG00000023236  1.426  7.067 245.0 3.133e-55 1.431e-51
## ENSMUSG00000050824  3.705  3.365 229.7 7.099e-52 2.882e-48
## ENSMUSG00000015852 -2.362  3.138 214.8 1.210e-48 4.049e-45
```


Get counts of significant genes (FDR-adjusted p-val < 0.01):

```r
tt.glm <- topTags(lrt, n = Inf)
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])
```

```
## [1] 600
```


Now grab the names of some highly significant genes:

```r
interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]
```

```
##                    SRX033480 SRX033488 SRX033481 SRX033489 SRX033482
## ENSMUSG00000020912     20.39    12.691     15.83    14.819    19.230
## ENSMUSG00000050141     14.14    10.153     10.67     6.264     9.946
## ENSMUSG00000035775     14.47    10.629     19.51    11.305    12.930
## ENSMUSG00000015484     26.97    33.949     32.76    36.206    31.497
## ENSMUSG00000024248     14.80    26.017     17.30    23.527    18.567
## ENSMUSG00000030532     27.96    24.906     21.35    25.360    18.567
## ENSMUSG00000054354     10.85     7.615     11.41     8.708     8.289
## ENSMUSG00000023236     65.45    78.684     72.50    72.413    76.587
##                    SRX033490 SRX033483 SRX033476 SRX033478 SRX033479
## ENSMUSG00000020912    17.330    26.430    13.359    23.621    14.310
## ENSMUSG00000050141     7.890     8.361     6.583    11.911    12.402
## ENSMUSG00000035775    12.540    12.136     6.389    13.527    11.210
## ENSMUSG00000015484    33.393    30.745    34.655    26.851    31.482
## ENSMUSG00000024248    25.784    18.879    24.781    14.334    14.310
## ENSMUSG00000030532    26.770    26.430    19.167    23.217    20.034
## ENSMUSG00000054354     7.608     7.012     3.291     4.442     5.963
## ENSMUSG00000023236    80.875    83.336    55.178    57.740    66.303
##                    SRX033472 SRX033473 SRX033474 SRX033475 SRX033491
## ENSMUSG00000020912    0.5037    0.4226    1.0805    0.2802    1.0995
## ENSMUSG00000050141    0.2518    0.0000    0.0000    0.2802    0.8246
## ENSMUSG00000035775    0.7555    0.0000    0.2701    0.2802    0.6872
## ENSMUSG00000015484    7.0516    4.8595    9.1841    7.8443    8.5209
## ENSMUSG00000024248    2.2666    1.2677    1.8908    3.3618    3.1610
## ENSMUSG00000030532   69.7605   70.5683   63.2082   63.0343   73.8023
## ENSMUSG00000054354    0.0000    0.0000    0.0000    0.0000    0.5497
## ENSMUSG00000023236  195.4301  183.6044  176.1186  156.3252  206.9762
##                    SRX033484 SRX033492 SRX033485 SRX033493 SRX033486
## ENSMUSG00000020912    0.0000    0.7027    0.8651    0.1362    0.0000
## ENSMUSG00000050141    0.0000    0.1405    0.0000    0.1362    0.2578
## ENSMUSG00000035775    0.6784    0.1405    0.5767    0.1362    1.8045
## ENSMUSG00000015484    9.4974    8.4319    8.9396    6.6759   10.3116
## ENSMUSG00000024248    1.5829    1.6864    1.7302    1.7712    2.5779
## ENSMUSG00000030532   70.0997   70.6873   62.5770   78.6123   59.2919
## ENSMUSG00000054354    0.0000    0.0000    0.0000    0.0000    0.0000
## ENSMUSG00000023236  177.9628  192.8090  207.0519  201.6399  204.6859
##                    SRX033494
## ENSMUSG00000020912    0.0000
## ENSMUSG00000050141    0.1477
## ENSMUSG00000035775    0.5907
## ENSMUSG00000015484    8.8604
## ENSMUSG00000024248    3.1012
## ENSMUSG00000030532   71.3265
## ENSMUSG00000054354    0.0000
## ENSMUSG00000023236  192.7144
```


Characterize the number of genes altered in the strains:

```r
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))
```

```
##    [,1] 
## -1   451
## 0  35660
## 1    425
```

There are 451 genes under-expressed in the DBA/2J mice, and 425 that are over-expressed. The remaining 35660 show no significant differences in expression.

Observe the fold-changes:

```r
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


DESeq analysis
--------------
Now try the DE analysis using the DESeq package.

Create the data set:

```r
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))
```

```
##                    SRX033480 SRX033488 SRX033481 SRX033489 SRX033482
## ENSMUSG00000000001       369       744       287       769       348
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         0         1         0         1         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         0         1         1         5         0
## ENSMUSG00000000049         0         1         0         1         0
##                    SRX033490 SRX033483 SRX033476 SRX033478 SRX033479
## ENSMUSG00000000001       803       433       469       585       321
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         1         0         7         6         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         4         0         0         0         0
## ENSMUSG00000000049         0         0         0         0         0
##                    SRX033472 SRX033473 SRX033474 SRX033475 SRX033491
## ENSMUSG00000000001       301       461       309       374       781
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         1         1         1         1         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         4         1         1         0         1
## ENSMUSG00000000049         0         0         0         0         0
##                    SRX033484 SRX033492 SRX033485 SRX033493 SRX033486
## ENSMUSG00000000001       555       820       294       758       419
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         2         1         1         4         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         2         1         1         1         1
## ENSMUSG00000000049         0         0         0         0         0
##                    SRX033494
## ENSMUSG00000000001       857
## ENSMUSG00000000003         0
## ENSMUSG00000000028         5
## ENSMUSG00000000031         0
## ENSMUSG00000000037         2
## ENSMUSG00000000049         0
```


Estimate size factors:

```r
deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)
```

```
## SRX033480 SRX033488 SRX033481 SRX033489 SRX033482 SRX033490 SRX033483 
##    0.6439    1.3454    0.5785    1.4295    0.6355    1.5240    0.7933 
## SRX033476 SRX033478 SRX033479 SRX033472 SRX033473 SRX033474 SRX033475 
##    1.1272    1.0772    0.8984    0.8886    1.0255    0.7987    0.7796 
## SRX033491 SRX033484 SRX033492 SRX033485 SRX033493 SRX033486 SRX033494 
##    1.6162    0.9882    1.5720    0.7558    1.5922    0.8264    1.4715
```

```r
deSeqDat <- estimateDispersions(deSeqDat)
```


Plot dispersion vs.counts:

```r
plotDispEsts(deSeqDat)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


Now fit the model:

```r
results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)
```

```
## 'data.frame':	36536 obs. of  8 variables:
##  $ id            : chr  "ENSMUSG00000000001" "ENSMUSG00000000003" "ENSMUSG00000000028" "ENSMUSG00000000031" ...
##  $ baseMean      : num  489.18 0 1.57 0 1.1 ...
##  $ baseMeanA     : num  509.685 0 1.657 0 0.859 ...
##  $ baseMeanB     : num  470.53 0 1.49 0 1.32 ...
##  $ foldChange    : num  0.923 NaN 0.898 NaN 1.537 ...
##  $ log2FoldChange: num  -0.115 NaN -0.156 NaN 0.62 ...
##  $ pval          : num  0.498 NA 0.829 NA 0.968 ...
##  $ padj          : num  1 NA 1 NA 1 ...
```

```r
plotMA(results)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 


Limma & voom analysis
---------------------
Use Voom and limma to perform DE analysis.

Prepare the data (figure out weights, etc.):

```r
norm.factor <- calcNormFactors(dat)
show(norm.factor)
```

```
##  [1] 0.9888 0.9783 1.0034 1.0041 0.9949 0.9850 0.9995 1.0059 0.9967 0.9873
## [11] 1.0250 0.9950 0.9923 1.0137 1.0067 1.0408 0.9961 1.0234 0.9816 0.9977
## [21] 0.9861
```

```r
dat.voomed <- voom(dat, design, plot = TRUE, lib.size = colSums(dat) * norm.factor)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

```r
dat.voomed
```

```
## An object of class "EList"
## $E
##                    SRX033480 SRX033488 SRX033481 SRX033489 SRX033482
## ENSMUSG00000000001     6.941     6.916     6.720     6.871     6.860
## ENSMUSG00000000003    -2.588    -3.625    -2.447    -3.716    -2.585
## ENSMUSG00000000028    -2.588    -2.040    -2.447    -2.132    -1.000
## ENSMUSG00000000031    -2.588    -3.625    -2.447    -3.716    -2.585
## ENSMUSG00000000037    -2.588    -2.040    -0.862    -0.257    -2.585
##                    SRX033490 SRX033483 SRX033476 SRX033478 SRX033479
## ENSMUSG00000000001    6.8447      6.87    6.4977    6.8899     6.279
## ENSMUSG00000000003   -3.8055     -2.89   -3.3773   -3.3036    -3.049
## ENSMUSG00000000028   -2.2205     -2.89    0.5296    0.3968    -1.465
## ENSMUSG00000000031   -3.8055     -2.89   -3.3773   -3.3036    -3.049
## ENSMUSG00000000037   -0.6355     -2.89   -3.3773   -3.3036    -3.049
##                    SRX033472 SRX033473 SRX033474 SRX033475 SRX033491
## ENSMUSG00000000001    6.2110     6.615     6.397     6.693     6.737
## ENSMUSG00000000003   -3.0250    -3.236    -2.877    -2.855    -3.873
## ENSMUSG00000000028   -1.4400    -1.651    -1.292    -1.270    -2.288
## ENSMUSG00000000031   -3.0250    -3.236    -2.877    -2.855    -3.873
## ENSMUSG00000000037    0.1449    -1.651    -1.292    -2.855    -2.288
##                    SRX033484 SRX033492 SRX033485 SRX033493 SRX033486
## ENSMUSG00000000001    6.9151     6.855     6.375     6.718     6.760
## ENSMUSG00000000003   -3.2025    -3.825    -2.827    -3.849    -2.952
## ENSMUSG00000000028   -0.8806    -2.240    -1.242    -0.679    -1.367
## ENSMUSG00000000031   -3.2025    -3.825    -2.827    -3.849    -2.952
## ENSMUSG00000000037   -0.8806    -2.240    -1.242    -2.264    -1.367
##                    SRX033494
## ENSMUSG00000000001     7.005
## ENSMUSG00000000003    -3.739
## ENSMUSG00000000028    -0.280
## ENSMUSG00000000031    -3.739
## ENSMUSG00000000037    -1.417
## 36531 more rows ...
## 
## $weights
##        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]  [,10]
## [1,] 15.907 19.293 15.443 19.597 15.899 19.891 16.895 18.482 18.240 17.413
## [2,]  1.632  1.547  1.632  1.531  1.632  1.517  1.632  1.595  1.612  1.632
## [3,]  1.469  1.383  1.488  1.378  1.470  1.374  1.436  1.397  1.402  1.421
## [4,]  1.632  1.547  1.632  1.531  1.632  1.517  1.632  1.595  1.612  1.632
## [5,]  1.558  1.421  1.585  1.414  1.558  1.407  1.507  1.445  1.453  1.484
##       [,11]  [,12]  [,13]  [,14]  [,15]  [,16]  [,17]  [,18]  [,19]  [,20]
## [1,] 16.988 17.670 16.503 16.432 19.760 17.564 19.603 16.339 19.682 16.748
## [2,]  1.632  1.632  1.632  1.632  1.531  1.632  1.539  1.632  1.535  1.632
## [3,]  1.403  1.390  1.415  1.416  1.354  1.392  1.358  1.419  1.356  1.409
## [4,]  1.632  1.632  1.632  1.632  1.531  1.632  1.539  1.632  1.535  1.632
## [5,]  1.423  1.406  1.437  1.439  1.370  1.409  1.373  1.442  1.371  1.429
##       [,21]
## [1,] 19.318
## [2,]  1.555
## [3,]  1.364
## [4,]  1.555
## [5,]  1.377
## 36531 more rows ...
## 
## $design
##   (Intercept) group2
## 1           1      0
## 2           1      0
## 3           1      0
## 4           1      0
## 5           1      0
## 16 more rows ...
## 
## $targets
## [1] 3006197 6166939 2726341 6572633 3000692
## 16 more rows ...
```


Now fit the model using our "voomed" data, and check out the top results:

```r
fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)
topTable(fit)
```

```
##                    X.Intercept.    group2 AveExpr      F   P.Value
## ENSMUSG00000025867       12.477 -0.013979  12.471 154234 2.584e-55
## ENSMUSG00000022892       12.268 -0.086479  12.223 138459 1.103e-54
## ENSMUSG00000037852       12.321 -0.237470  12.194 134826 1.576e-54
## ENSMUSG00000042700       10.716 -0.048029  10.692 133973 1.717e-54
## ENSMUSG00000029461       10.300  0.020715  10.311 121795 6.182e-54
## ENSMUSG00000020658        9.610 -0.019628   9.601 119935 7.604e-54
## ENSMUSG00000060261        9.469 -0.015743   9.461 117841 9.635e-54
## ENSMUSG00000032549       11.904  0.003545  11.905 117324 1.022e-53
## ENSMUSG00000024462       10.227 -0.138929  10.153 116767 1.090e-53
## ENSMUSG00000030102       12.085 -0.026149  12.073 112155 1.874e-53
##                    adj.P.Val
## ENSMUSG00000025867 9.442e-51
## ENSMUSG00000022892 1.568e-50
## ENSMUSG00000037852 1.568e-50
## ENSMUSG00000042700 1.568e-50
## ENSMUSG00000029461 4.424e-50
## ENSMUSG00000020658 4.424e-50
## ENSMUSG00000060261 4.424e-50
## ENSMUSG00000032549 4.424e-50
## ENSMUSG00000024462 4.424e-50
## ENSMUSG00000030102 6.417e-50
```


Take-home
---------
Choose a p-value cutoff, and identify DE genes from each of the 3 methods:

```r
cutoff <- 1e-05
edger.hits <- rownames(tt.glm$table[tt.glm$table$FDR < cutoff, ])
deseq.hits <- results[which(results$padj < cutoff), "id"]
voom.hits <- rownames(topTable(fit, number = Inf, coef = "group2", p.value = cutoff))
```


Plot the overlaps using a venn diagram:

```r
library(VennDiagram)
```

```
## Loading required package: grid
```

```r
de.genes <- list(edgeR = edger.hits, DESeq = deseq.hits, Voom = voom.hits)
venn.plot <- venn.diagram(de.genes, filename = NULL, fill = c("red", "blue", 
    "green"), cat.fontface = "bold", cat.cex = 1.5, cex = 2, fontface = "bold")
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 


Plot the counts for a few of the shared genes:

```r
common.genes <- intersect(edger.hits, intersect(deseq.hits, voom.hits))
set.seed(540)
my.genes <- common.genes[sample(c(1:length(common.genes)), size = 3)]
my.dat <- dat[my.genes, ]
my.dat <- data.frame(counts = as.vector(t(my.dat)), gene = rep(rownames(my.dat), 
    each = ncol(my.dat)))
my.dat <- suppressWarnings(data.frame(des, my.dat))
bwplot(counts ~ strain | gene, my.dat, panel = panel.violin, scales = list(y = list(relation = "free")), 
    col = "darkgrey", border = "darkgrey")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 

