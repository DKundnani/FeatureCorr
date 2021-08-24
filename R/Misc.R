library(preprocessCore)
ztransfun<-function(col){
  (col-mean(as.numeric(unlist(col))))/sd(as.numeric(unlist(col)))
}
#' Data Transformation and filtration
#'
#' Transformation of data using one of the four mentioned methods and filtering based on median value
#' @param df numerical dataframe with rows having series of values for a single feature
#' @param transformation type of transformation method ''Log2' or 'Z-score' or 'quantile' or 'NA' (no) transformaiton
#' @param margin 1 indicates rows, 2 indicates columns
#' 
#' @examples
#' transdf <- data_transform(df=GTEXv7[,3:10],transformation='log2', marg=2);
#' @docType data
#'
#' @usage data(GTEXv7)
#' 
#' @export
data_transform <- function(df,transformation='log2', marg=2){ 
  
  if (marg==1){ df=as.data.frame(t(df))}
  
  if (transformation == 'log2') {
    newdf=log2(df+1)
  } else if (transformation == 'Z-score') {
    newdf<-apply(df, 2, ztransfun)
  } else if (transformation == 'quantile') {
    newdf<-normalize.quantiles(as.matrix(df))
  } else if (transformation == 'NA') {
    newdf<-normalize.quantiles(as.matrix(df))
  } else {stop("There is no such transformation available")}
  
  colnames(newdf)<-colnames(df)
  
  OldMedianValue<-apply(df, 1, median, na.rm=F)
  NewMedianValue<-apply(newdf, 1, median, na.rm=F)
  Median=data.frame("MedianValues of Original DAtaframe" = OldMedianValue, 
                    "MedianValues of Transformed DAtaframe" = NewMedianValue)
  
  plot.new()
  par(mfrow=c(2, 2), mar=c(15,5,1,1))
  Om=hist(OldMedianValue, breaks = 100, xlab = "Median Value", col="chocolate", main ="Distribution of Median Value in Original Dataframe")
  On=hist(NewMedianValue, breaks = 100, xlab = "Median Value", col="chocolate", main ="Distribution of Median Value in Transformed Dataframe")
  OG.box<-boxplot(df, main ="Boxplots of Original Data", xaxt = "n")
  axis(side = 1,las = 2, labels = FALSE)
  text(x = 1:ncol(df),y = par("usr")[3] - 0.45,labels = colnames(df), srt = 90, adj=1.1, xpd = NA)
  new.box<-boxplot(newdf, main ="Boxplots of Transformed Data", xaxt = "n")
  axis(side = 1,las = 2, labels = FALSE)
  text(x = 1:ncol(df),y = par("usr")[3] - 0.45,labels = colnames(df), srt = 90, adj=1.1, xpd = NA)
  dist.plot<- recordPlot()
  
  if (marg==1){ newdf=as.data.frame(t(newdf))}
  
  return(list(newdf,Median,dist.plot))
  
}





