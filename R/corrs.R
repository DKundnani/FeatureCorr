#Corr
library(dplyr)
library(rlang)
library(Hmisc)
library(corrplot)
library(ggpubr)
library(MASS)
library(ggExtra)
library(GeneNet)
#' Prime Feature Correlation
#'
#' Get Correlation Coefficient and its distribution of a single prime feature with all the other features
#' @param df A dataframe with echo row corresponding to a feature
#' @param featurelist listing of feature names or ids to measure correlation between
#' @param primefeature one feature of interest to measure correlation against all other features
#' 
#' @examples
#' primefcorr <- primefeature_corr(df=GTEXv7[,3:ncol(GTEXv7)],featurelist=GTEXv7[,2] ,primefeature="RNASEH2A");
#' @docType data
#'
#' @usage data(GTEXv7)
#' 
#' @export
primefeature_corr <- function(df,featurelist,primefeature, quantiles=c(0,0.01,.01,.05,0.10,.25,.50,.75,.90,.95,.99,.999,1)){
  if (length(featurelist)!=nrow(df) && length(primefeature)>1) { stop("Error: Make sure the length of feature list is same as the number of rows in the dataframe and there is only one prime feature.n", call. = FALSE) }
  c=ncol(df)
  r=nrow(df)
  spearman_corr<-rep(NA, r)
  spearman_pval<-rep(NA, r)
  pearson_corr<-rep(NA, r)
  pearson_pval<-rep(NA, r)
  
  x=match(primefeature, featurelist)
  if (is.na(x)) {  stop("Primefeature not found in feature list.n", call. = FALSE)}
  
  writeLines("Calculating Correlation Coefficients")
  for (i in 1:r){
    s<-cor.test(as.numeric(df[i,]), as.numeric(df[x,]), method='spearman')
    spearman_corr[i]<-s$estimate
    spearman_pval[i]<-s$p.value
    
    p<-cor.test(as.numeric(df[i,]), as.numeric(df[x,]), method='pearson')
    pearson_corr[i]<-p$estimate
    pearson_pval<-p$p.value
  }
  
  df$Features<-featurelist
  df$Spearman_correlation_Coefficient<-spearman_corr
  df$Spearman_correlation_pvalue<-spearman_pval
  df$Pearson_correlation_Coefficient<-pearson_corr
  df$Pearson_correlation_pvalue<-pearson_pval
  correlations<-df[,(c+1):ncol(df)]
  
  writeLines("Calculating Distribution Statistics")
  quantiles_Spearman.df<-quantile(df$Spearman_correlation_Coefficient, quantiles)
  quantiles_Pearson.df<-quantile(df$Pearson_correlation_Coefficient, quantiles)
  
  quant=data.frame("Speaman Correlation Coefficient Quantiles" = quantiles_Spearman.df, 
                   "Pearson Correlation Coefficient Quantiles" = quantiles_Pearson.df)
  
  writeLines("Plotting Correlation _Coefficient Distribution")
  
  par(mfrow=c(2, 1), mar=c(7,5,1,1))
  s=hist(df$Spearman_correlation_Coefficient, breaks = 20, xlab = "Spearman Correlation Coefficient", col="steelblue3", main ="Distribution of Spearman Correlation Coefficient")
  p=hist(df$Pearson_correlation_Coefficient, breaks = 20, xlab = "Pearson Correlation Coefficient", col="seagreen4", main ="Distribution of Pearson Correlation Coefficient")
  freq.plot<- recordPlot()
  
  finaloutput<-list(df,correlations, quant,s,p,freq.plot)
  
  names(finaloutput) <- c("DatawithCorr", "Correlations","Quantiles","Spearmann_Corr.freq","Pearson.freq",  "Frequency_Dist")
  return(finaloutput)
}

#' Pairwise Correlation
#'
#' Get Pair wise correlations for a given list of features 
#' @param df numerical dataframe with rows having series of values for a single feature
#' @param featurelist list of features in the sequence as it occurs in the Dataframe
#' @param featuregroup list of TRUE/FALSE values for features to be observed vs rest of the features
#' @param clusterorder Character, the ordering method of the correlation matrix.'original' for original order (default),'AOE' for the angular order of the eigenvectors, 'FPC' for the first principal component order, 'hclust' for the hierarchical clustering order, 'alphabet' for alphabetical order.
#' @param clustermethod Character, the ordering method of the correlation matrix.'original' for original order (default),'AOE' for the angular order of the eigenvectors, 'FPC' for the first principal component order, 'hclust' for the hierarchical clustering order, 'alphabet' for alphabetical order.
#' 
#' @examples
#' pairwisedf<-t(TCGA40[,3:ncol(TCGA40)]) 
#' featgroup<-grepl( "RNASE",rownames(pairwisedf)) #optional
#' pairwisecorr <- pairwise_corr(df=pairwisedf,featurelist=rownames(pairwisedf), featuregroup=featgroup, clusterorder="hclust");
#' @docType data
#'
#' @usage data(TCGA40)
#' 
#' @export
pairwise_corr <- function(df,featurelist, featuregroup='NA', clusterorder="hclust"){
  if (length(featurelist)!=nrow(df)) { stop("Error: Make sure the length of feature list is same as the number of rows in the dataframe and there is only one prime feature.n", call. = FALSE) }
  rownames(df)<-featurelist
  corr_pair<-rcorr(t(as.matrix(df)))
  par(mfrow=c(1, 1), mar=c(1,4,4,8))
  if (featuregroup == 'NA') {
    corrplot(corr_pair$r,type = "full", order =  "hclust",addrect = 2,method = "square",tl.col = "black",tl.srt = 45, cl.cex = 1.3, pch.cex = 0.5)
    out=list(pairwise_corr$r,pairwise_corr$P)
  } else {
    corrplot(corr_pair$r[featuregroup,],method = "square",tl.col = "black",tl.srt = 90, cl.cex = 1.3, pch.cex = 0.5, is.corr = FALSE, cl.pos = 'b', cl.ratio = 0.6,)
    out=list(pairwise_corr$r,pairwise_corr$P)
  }
  return(out)
  
}

#' Pair Scatter
#'
#' Scatter between two chosen features
#' @param df numerical dataframe with rows having series of values for a single feature
#' @param featurelist list of features in the sequence as it occurs in the Dataframe
#' @param feature1,feature2 Features to be compared
#' 
#' @docType data
#' 
#' @usage data(TCGA40)
#' 
#' @examples
#' featgroup<-grepl( "RNASE",rownames(TCGA40)) #optional
#' pscatter <- pair_scatter(df=TCGA40,featurelist=rownames(TCGA40), feature1="RNASEH2A", feature2="RNASEH2C");
#' @export
pair_scatter <- function(df,featurelist, feature1, feature2){
  if (length(featurelist)!=nrow(df)) { stop("Error: Make sure the length of feature list is same as the number of rows in the dataframe and there is only one prime feature.n", call. = FALSE) }
  
  data=as.data.frame(t(df))
  colnames(data)<-featurelist
  #scatter plot with Pearson's R
  par(mfrow = c(1, 1))
  p<-ggscatter(data, x = feature1 , y = feature2, color = "#2c7fb8", margin.params = list(fill = "lightgray"), add = "reg.line", add.params = list(color='red', fille='lightgray'),conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method='pearson',label.x.npc = 0,label.y.npc = 1, size=10)) +theme_bw()
  ggMarginal(p, size = 2, type = "histogram", col = "blue", fill = "#2c7fb8")
  p <- recordPlot()
  
  plot.new()
  kern <- kde2d(as.numeric(unlist(data[feature1])),as.numeric(unlist(data[feature2])))
  smoothScatter(as.numeric(unlist(data[feature1])),as.numeric(unlist(data[feature2])), xlab=feature1, ylab =feature2, nrpoints = 0)
  contour(kern, drawlabels = FALSE, nlevels = 6,
          col = rev(heat.colors(6)), add = TRUE, lwd = 3)
  q <- recordPlot()
  ## clean up device
  q
  return(list(p,q))
}  

















