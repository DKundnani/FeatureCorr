ztransfun<-function(col){
  (col-mean(as.numeric(unlist(col))))/sd(as.numeric(unlist(col)))
}
#' Data Transformation and filtration
#' Transformation of data using one of the four mentioned methods and filtering based on median value
#'
#'
#' @param df numerical dataframe with rows having series of values for a single feature
#' @param transformation type of transformation method ''Log2' or 'Z-score' or 'quantile' or 'NA' (no) transformaiton
#' @param featurelist listing of feature names or ids to measure correlation between
#' @param medianthres (optional)features have median below this numerical value are filtered or removed
#'
#' @usage data_transform(df,transformation,featurelist, medianthres)
#'
#' @examples
#' transdf<- data_transform(df=GTEX[-1],transformation='log2',
#'                          featurelist=GTEX$Description)
#' transdf<- data_transform(df=GTEX[-1],transformation='log2',
#'                          featurelist=GTEX$Description, medianthres=1)
#'
#' @export

data_transform <- function(df,transformation='log2',featurelist, medianthres='NA'){
  if (transformation == 'log2') {
    newdf=log2(df+1)
  } else if (transformation == 'Z-score') {
    newdf<-apply(df, 2, ztransfun)
  } else if (transformation == 'quantile') {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(c("preprocessCore"))
    newdf<-normalize.quantiles(as.matrix(df))
  } else if (transformation == 'NA') {
    newdf<-normalize.quantiles(as.matrix(df))
  } else {stop("There is no such transformation available")}

  colnames(newdf)<-colnames(df)
  rownames(newdf)<-featurelist

  OldMedianValue<-apply(df, 1, median, na.rm=F)
  NewMedianValue<-apply(newdf, 1, median, na.rm=F)
  Median=data.frame("MedianValues of Original DAtaframe" = OldMedianValue,
                    "MedianValues of Transformed DAtaframe" = NewMedianValue)

  if (medianthres!='NA'){
    newdf=newdf[NewMedianValue>=medianthres,]
    NewMedianValue<-NewMedianValue[NewMedianValue>=medianthres]
  }

  if (ncol(df)<15) {
    par(mfrow=c(2, 2), mar=c(15,5,2,2))
  } else {
    par(mfrow=c(1, 2), mar=c(5,5,2,2))
  }
  Om=hist(OldMedianValue, breaks = 200, xlab = "Median Value", col="chocolate",lty="blank", main ="Original Data")
  On=hist(NewMedianValue, breaks = 200, xlab = "Median Value", col="chocolate",lty="blank", main ="Transformed Data")

  if (ncol(df)<15) {
    OG.box<-boxplot(df, main ="Boxplots of Original Data", xaxt = "n")
    axis(side = 1,las = 2, labels = FALSE)
    text(x = 1:ncol(df),y = par("usr")[3] - 0.45,labels = colnames(df), srt = 90, adj=1.1, xpd = NA)

    if (medianthres!='NA') {
      new.box<-boxplot(newdf[ , -which(names(newdf) %in% c("feature"))], main ="Boxplots of Transformed Data", xaxt = "n")
      axis(side = 1,las = 2, labels = FALSE)
      text(x = 1:ncol(df),y = par("usr")[3] - 0.45,labels = colnames(df), srt = 90, adj=1.1, xpd = NA)
    } else {
      new.box<-boxplot(newdf, main ="Boxplots of Transformed Data", xaxt = "n")
      axis(side = 1,las = 2, labels = FALSE)
      text(x = 1:ncol(df),y = par("usr")[3] - 0.45,labels = colnames(df), srt = 90, adj=1.1, xpd = NA)
    }

  }

  dist.plot<- recordPlot()

  return(list(newdf,Median,dist.plot))

}





