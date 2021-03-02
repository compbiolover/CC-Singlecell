#Name: mad_calculator.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently build median absolute
#deviation (MAD) metric gene lists. 

#The median absolute deviation (MAD) metric ----
mad_calculator <- function(denoised.sc=sc){
  mads <- apply(denoised.sc,1,mad)
  index <- order(mads, decreasing=TRUE)
  mad.ranking<- mads[index]
  mad.ranking<-abs(mad.ranking)/sum(abs(mad.ranking))
  
  #Returning the finished MAD list----
  return(mad.ranking)
}
