
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FeatureCorr

<!-- badges: start -->

<!-- badges: end -->

An R package to study feature correlations aided with data
transformation for Next Generation sequencing and microarray data

## Installation

You can install the released version of FeatureCorr from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("FeatureCorr")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DKundnani/FeatureCorr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(FeatureCorr) #Loading the library
```

Data Transformation using FeatureCorr::data\_transform function

``` r
#Data Transformation
data("GTEX") #Load GTEX data names GTEXv7
transdf<- data_transform(df=GTEXv7[-1],transformation='log2', featurelist=GTEXv7$Description, medianthres=1)
```

<img src="man/figures/README-preprocessing-1.png" width="100%" />

Visualizing Single Pair Correlation using using
FeatureCorr::pair\_scatter function

``` r
data("TCGA40") #Load TCGAdata for 40 genes
data("logTCGA40") #Load log transformed TCGAdata for 40 genes
#Single Pair Scatter on Original RNA-seq Dataset vs log2 transformed RNA-seq Dataset
pscatter <- pair_scatter(df=TCGA40,featurelist=rownames(TCGA40),feature1="RNASEH2A",feature2="PCNA", corrmeth='pearson')
#> `geom_smooth()` using formula 'y ~ x'
#> `geom_smooth()` using formula 'y ~ x'
```

<img src="man/figures/README-checking-pairwise-scatter-1.png" width="100%" />

``` r

pscatter <- pair_scatter(df=logTCGA40,featurelist=rownames(logTCGA40),feature1="RNASEH2A", feature2="PCNA", corrmeth='pearson')
#> `geom_smooth()` using formula 'y ~ x'
#> `geom_smooth()` using formula 'y ~ x'
```

<img src="man/figures/README-checking-pairwise-scatter-2.png" width="100%" />

Data Transformation using FeatureCorr::primefeature\_corr function

``` r
data("GTEX") #Load GTEX data names GTEXv7
#Data transformation and Prime Feature(RNASEH2A gene) Correlation in Gene-Tissue EXpression Dataset
transdf<- data_transform(df=GTEXv7[-1],transformation='log2', featurelist=GTEXv7$Description,
                          medianthres=0.5)
```

<img src="man/figures/README-prime-feature-correlation-workflow-1.png" width="100%" />

``` r
inputdf<-transdf[[1]] #First item in the list returned is the transformed Data
primefcorr <- primefeature_corr(df=inputdf,featurelist=rownames(inputdf) ,primefeature="RNASEH2A")
#> Calculating Correlation Coefficients
#> Calculating Distribution Statistics
#> Plotting Correlation Coefficient Distribution
```

<img src="man/figures/README-prime-feature-correlation-workflow-2.png" width="100%" />

Visualizing Multiple Pair-wise Correlation using
FeatureCorr::pairwise\_corr function

``` r
data("TCGA40") #Load TCGAdata for 40 genes
data("logTCGA40") #Load log transformed TCGAdata for 40 genes
#Multiple pairwise correlation for 40 genes in ~10K samples in The Cancer Genome Atlas-Pan Cancer dataset
pairwisedf<-TCGA40[,3:ncol(TCGA40)] #Data provided with the Package
corr_pair<-pairwise_corr(df=logTCGA40,featurelist=rownames(logTCGA40), visorder="hclust")
```

<img src="man/figures/README-multiple-pariwise-correlation-1.png" width="100%" />

``` r

# Separating feature of interest for visualization 
# In This technique you can group features of differen types. For example if I have methylation and transcription data in the same dataframe for same set of feature identfiers
featgroup<-grepl( "RNASE",rownames(logTCGA40)) #optional, a set of features to separated
corr_pair<-pairwisecorr <- pairwise_corr(df=pairwisedf,featurelist=rownames(pairwisedf),featuregroup=featgroup)
```

<img src="man/figures/README-multiple-pariwise-correlation-2.png" width="100%" />
