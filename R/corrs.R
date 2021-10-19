#Corr

#' Prime Feature Correlation
#' Get Correlation Coefficient and its distribution of a single prime feature with all the other features
#'
#' @import dplyr
#' @import rlang
#' @importFrom Hmisc rcorr
#' @import corrplot
#' @import ggpubr
#' @importFrom MASS kde2d
#' @import ggExtra
#' @import GeneNet
#'
#'
#' @param df A dataframe with echo row corresponding to a feature
#' @param featurelist listing of feature names or ids to measure correlation between
#' @param primefeature one feature of interest to measure correlation against all other features
#' @param corrmeth Correlation method used. 'pearson' or 'spearman' (default='pearson')
#' @param quantiles list of corresponding probabilities for quantiles to be obtained
#' (default = c(0,0.01,.01,.05,0.10,.25,.50,.75,.90,.95,.99,.999,1) )
#'
#' @examples
#' transdf<- data_transform(df=GTEX[-1],transformation='log2',
#'                          featurelist=GTEX$Description, medianthres=1)
#' inputdf<-transdf[[1]] #First item in the list returned is the transformed Data
#' primefcorr <- primefeature_corr(df=inputdf,featurelist=rownames(inputdf),primefeature="RNASEH2A");
#'
#' @usage primefeature_corr(df,featurelist,primefeature,corrmeth,quantiles)
#'
#' @export
primefeature_corr <- function(df,featurelist,primefeature, corrmeth='pearson', quantiles=c(0,0.01,.01,.05,0.10,.25,.50,.75,.90,.95,.99,.999,1)){
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
#' @param visorder  The ordering method of the correlation matrix.'original' for original order ,'AOE' for the angular order of the eigenvectors, 'FPC' for the first principal component order, 'hclust'(default) for the hierarchical clustering order, 'alphabet' for alphabetical order.
#' @param corrmeth Correlation method used. 'pearson' or 'spearman' (default='pearson')
#' @param clustno number of clusters expected to be drawn on the plot (default = NULL)
#' @examples
#' corr_pair<-pairwise_corr(df=logTCGA40,featurelist=rownames(logTCGA40), visorder="hclust", clustno=2)
#' featgroup<-grepl( "RNASE",rownames(logTCGA40)) #optional, a set of TRUE/FALSE of length featurelist
#' corr_pair<-pairwisecorr <- pairwise_corr(df=logTCGA40,featurelist=rownames(logTCGA40),
#'                                          featuregroup=featgroup)
#'
#' @usage pairwise_corr(df,featurelist, featuregroup,visorder,corrmeth, clustno=NULL)
#'
#' @export
pairwise_corr <- function(df,featurelist, featuregroup='NA', visorder="hclust", corrmeth='pearson', clustno = NULL){
  if (length(featurelist)!=nrow(df)) { stop("Error: Make sure the length of feature list is same as the number of rows in the dataframe", call. = FALSE) }
  rownames(df)<-featurelist
  corr_pair<-rcorr(t(as.matrix(df)), type=corrmeth)
  par(mfrow=c(1, 1), mar=c(1,4,4,8))

  if (length(featuregroup) == 1) {
    corrplot(corr_pair$r,type = "full", order =  "hclust",addrect = clustno ,method = "square",tl.col = "black",tl.srt = 60, cl.cex = 1.3, pch.cex = 0.5, tl.cex = 0.7)
    out=list(corr_pair$r,corr_pair$P)

  } else {
    corrplot(corr_pair$r[featuregroup,!featuregroup],method = "square",tl.col = "black",tl.srt = 90, cl.cex = 1.3, pch.cex = 0.5, is.corr = FALSE, cl.pos = 'b', cl.ratio = 0.6)
    out=list(corr_pair$r,corr_pair$P)

  }
  return(out)

}

#' Pair Scatter
#' Scatter between two chosen features.
#'
#' @import ggpubr
#' @import ggplot2
#'
#'
#' @param df numerical dataframe with rows having series of values for a single feature
#' @param featurelist list of features in the sequence as it occurs in the Dataframe
#' @param feature1,feature2 Features to be compared
#' @param corrmeth Correlation method used. 'pearson' or 'spearman' (default='pearson')
#' @param col Color for scatter plot, default="#1f78b4"
#'
#' @usage pair_scatter(df,featurelist, feature1, feature2, corrmeth)
#'
#' @examples
#' pscatter <- pair_scatter(df=TCGA40,featurelist=rownames(TCGA40),
#'                          feature1="RNASEH2A",feature2="PCNA", corrmeth='pearson', col= "#1f78b4")
#' pscatter <- pair_scatter(df=logTCGA40,featurelist=rownames(logTCGA40),
#'                          feature1="RNASEH2A", feature2="PCNA", corrmeth='pearson', col= "#b45b1f")
#'
#' @export
pair_scatter <- function(df,featurelist, feature1, feature2,corrmeth='pearson', col= "#1f78b4"){
  if (length(featurelist)!=nrow(df)) { stop("Error: Make sure the length of feature list is same as the number of rows in the dataframe", call. = FALSE) }

  data=as.data.frame(t(df))
  colnames(data)<-featurelist
  #scatter plot
  #par(mfrow = c(1, 1))
  o <- ggscatter(data, x = feature1 , y = feature2,  size = 0.7, color = col, margin.params = list(fill = "lightgray"), add = "reg.line", add.params = list(color='red', fille='lightgray'),conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method=corrmeth,label.x.npc = 0,label.y.npc = 1, size=10)) +theme_bw()
  p <-ggMarginal(o, size = 2, type = "histogram", col = "blue", fill = "#2c7fb8")

  plot.new()
  kern <- kde2d(as.numeric(unlist(data[feature1])),as.numeric(unlist(data[feature2])))
  smoothScatter(as.numeric(unlist(data[feature1])),as.numeric(unlist(data[feature2])), xlab=feature1, ylab =feature2, nrpoints = 0)
  contour(kern, drawlabels = FALSE, nlevels = 6,
          col = rev(heat.colors(6)), add = TRUE, lwd = 3)
  q <- recordPlot()

  ## clean up device
  plot.new()

  return(list(p,q))
}



