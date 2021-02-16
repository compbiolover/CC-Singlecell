databin<- function(dat){
  databinary<- dat
  databinary[databinary>0] <- 1
  databinary[databinary<0] <- 0
  return(databinary)
}

dataclean <- function(dat){
  message("Remove genes that all are zeros...")
  dat <- dat[rowSums(dat)!=0,]
  message("done")
  return(dat)
}

idfyImpgene <- function(data,databinary,condition,thres=0.05){

  
  beta_glm <- bayeswald(databinary,condition)
  genes_fit <- rownames(beta_glm)[beta_glm[,"chi-test"]< thres]
  
 # genes_fit<- genes_fit[! (genes_fit %in% rownames(databinary[rowSums(databinary)<2 ,]))]
  return(genes_fit)
}

idfyUsegene <- function(data,databinary,condition,ratio=0.8){

  genes_use<-rownames(databinary[rowSums(databinary)>=(ncol(data)*ratio), ])
  #genes_use<- genes_use[!genes_use %in% genes_fit]
  return(genes_use)
}

FunImpute <- function(
  object,
  genes_use = NULL,
  genes_fit = NULL,
  s.use = 20,
  do.print = FALSE,
  gram = TRUE,
  dcorgene = NULL
) {
  
  ifelse(is.null(genes_use),genes_use<- (rownames(object)),genes_use)
  ifelse(is.null(genes_fit),genes_fit<-rownames(object),genes_fit)
  ifelse(is.null(dcorgene),dcorgene<- min(100,floor(length(genes_use)/2) ) ,dcorgene)
  
  tsip<-object[genes_fit,]
  tsg<- object[genes_use,]
  
  corgene<- sapply(
    data.frame(t(tsip)), FUN = function(x,x2){
      
      genecorre <- apply( x2[which( x!=0),], 2, FUN = function(x,b){
        ifelse(length(unique(x)==1) ,x[1]<-x[1]+1e-5 ,x )
        ifelse(length(unique(b)==1) ,b[1]<-b[1]+1e-5 ,b )
        cor(x,b)
      },  x[which( x!=0)])
      
      genes_use[order(-genecorre)[1:(dcorgene)]]
      
    },t(tsg)
  )
  colnames(corgene)<-genes_fit
  corgene<- rbind(colnames(corgene),corgene)
 
  lasso_fits<- sapply(data.frame(corgene), FUN = function(x,tsip,tsg){
    
    x<-as.character(x)
    cell_loc <- which(tsip[x[1],]==0)
    
    lasso_input <- tsg[x[-1] , which(tsip[ x[1], ]!=0) ]
    lasso_fituse <- tsip[x[1] , which(tsip[ x[1], ]!=0) ]
    gene_loc <- match(rownames(lasso_input),rownames(tsg))

    lasso_fit<- sapply(data.frame(tsg[ ,cell_loc] ),FUN = function(x,lasso_input,lasso_fituse,gene_loc){
      lasso_out <- x[gene_loc]
      ##if genes x[-1] in tsg are unique_length==1,then fit<-0
      
      ifelse((length(unique(lasso_out))==1)|( length(unique(lasso_input)) == length( lasso_fituse ) ) , 0 , (lasso.fxn(
        lasso.input = lasso_input,
        lasso.output = lasso_out,
        lasso.fituse = lasso_fituse,
        s.use = s.use,
        do.print = do.print,
        gram = gram
      ))  )
      
    }, lasso_input,lasso_fituse,gene_loc)
    
  }, tsip,tsg)
  names(lasso_fits)<- colnames(corgene)
  
  imputedobj<- object
  for (i in 1:length(lasso_fits)){
    imputedobj[names(lasso_fits)[i], which(imputedobj[names(lasso_fits)[i] ,] ==0)] <- lasso_fits[[i]]
  }
  imputedobj[imputedobj<0] <- 0
  return(list(alldat =imputedobj,fitHur=lasso_fits))
}

lasso.fxn <- function(
  lasso.input,
  lasso.output,
  lasso.fituse,
  s.use = 20,
  gene.name = NULL,
  do.print = FALSE,
  gram = TRUE
) {
  lasso.model <- lars(
    x = lasso.input,
    y = as.numeric(x = lasso.output),
    type = "lasso",
    max.steps = s.use * 2,
    use.Gram = gram
  )
  #lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=min(s.use,max(lasso.model$df)))$fit
  lasso.coef <- coef.lars(
    object = lasso.model
  )[length(lasso.model$df),]
  
  lasso.fits<- lasso.fituse %*% (as.matrix(lasso.coef))
  #t(as.matrix(lasso.fituse[2,])) %*% ((lasso.coef))
  
  return(lasso.fits)
}

