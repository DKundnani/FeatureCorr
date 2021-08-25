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
#' transdf<- data_transform(df=GTEXv7[,3:ncol(GTEXv7)],transformation='log2', featurelist=GTEXv7[,2], medianthres=0.5)
#' inputdf<-transdf[[1]] #First item in the list returned is the transformed Data
#' primefcorr <- primefeature_corr(df=inputdf[,-ncol(inputdf)],featurelist=inputdf$feature ,primefeature="RNASEH2A");
#' @docType data
#'
#' @usage data(GTEXv7)
#'
#' @export
#'
primefeature_corr <- function(df,featurelist,primefeature,corrmeth='pearson', quantiles=c(0,0.01,.01,.05,0.10,.25,.50,.75,.90,.95,.99,.999,1)){
  if (length(featurelist)!=nrow(df) && length(primefeature)>1) { stop("Error: Make sure the length of feature list is same as the number of rows in the dataframe and there is only one prime feature.n", call. = FALSE) }
  c=ncol(df)
  r=nrow(df)
  corr<-rep(NA, r)
  pval<-rep(NA, r)


  x=match(primefeature, featurelist)
  if (is.na(x)) {  stop("Primefeature not found in feature list.n", call. = FALSE)}

  writeLines("Calculating Correlation Coefficients")
  for (i in 1:r){

    p<-cor.test(as.numeric(df[i,]), as.numeric(df[x,]), method=corrmeth)
    corr[i]<-p$estimate
    pval<-p$p.value
  }

  df$Features<-featurelist
  df$Correlation_Coefficient<-corr
  df$Correlation_pvalue<-pval

  correlations<-df[,(c+1):ncol(df)]

  writeLines("Calculating Distribution Statistics")

  quantiles<-quantile(df$Correlation_Coefficient, quantiles)


  writeLines("Plotting Correlation Coefficient Distribution")

  par(mfrow=c(1, 1), mar=c(7,5,1,1))
  s=hist(df$Correlation_Coefficient, breaks = 200, xlab = "Correlation Coefficient", col="steelblue3", main ="Distribution of Correlation Coefficient")

  freq.plot<- recordPlot()

  finaloutput<-list(df,correlations, quantiles,s,freq.plot)

  names(finaloutput) <- c("DatawithCorr", "Correlations","Quantiles","Corr.freq",  "Frequency_Dist")
  return(finaloutput)
}

#' Pairwise Correlation
#'
#' Get Pair wise correlations for a given list of features
#' @param df numerical dataframe with rows having series of values for a single feature
#' @param featurelist list of features in the sequence as it occurs in the Dataframe
#' @param featuregroup (optional)list of TRUE/FALSE values for features to be observed vs rest of the features(ordering/clustering is disabled)
#' @param clusterorder Character, the ordering method of the correlation matrix.'original' for original order (default),'AOE' for the angular order of the eigenvectors, 'FPC' for the first principal component order, 'hclust' for the hierarchical clustering order, 'alphabet' for alphabetical order.
#' @param clustermethod Character, the ordering method of the correlation matrix.'original' for original order (default),'AOE' for the angular order of the eigenvectors, 'FPC' for the first principal component order, 'hclust' for the hierarchical clustering order, 'alphabet' for alphabetical order.
#'
#' @examples
#' pairwisedf<-TCGA40[,3:ncol(TCGA40)]
#' corr_pair<-pairwise_corr(df=pairwisedf,featurelist=rownames(pairwisedf), clusterorder="hclust")
#' featgroup<-grepl( "RNASE",rownames(pairwisedf)) #optional, a set of features to separated
#' corr_pair<-pairwisecorr <- pairwise_corr(df=pairwisedf,featurelist=rownames(pairwisedf), featuregroup=featgroup);
#' @docType data
#'
#' @usage data(TCGA40)
#'
#' @export
pairwise_corr <- function(df,featurelist, featuregroup='NA', clusterorder="hclust"){
  if (length(featurelist)!=nrow(df)) { stop("Error: Make sure the length of feature list is same as the number of rows in the dataframe", call. = FALSE) }
  rownames(df)<-featurelist
  corr_pair<-rcorr(t(as.matrix(df)))
  par(mfrow=c(1, 1), mar=c(1,4,4,8))
  if (featuregroup == 'NA') {
    corrplot(corr_pair$r,type = "full", order =  "hclust",addrect = 2,method = "square",tl.col = "black",tl.srt = 45, cl.cex = 1.3, pch.cex = 0.5)
    out=list(corr_pair$r,corr_pair$P)
  } else {
    corrplot(corr_pair$r[featuregroup,!featuregroup],method = "square",tl.col = "black",tl.srt = 90, cl.cex = 1.3, pch.cex = 0.5, is.corr = FALSE, cl.pos = 'b', cl.ratio = 0.6,)
    out=list(corr_pair$r,corr_pair$P)
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

















