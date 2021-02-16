EMD_nonzero <- function(data.df, outcomes, emd, emd.perm) { 
  
  structure(list("data"=data.df, "outcomes"=outcomes,
                 "emd"=emd, "emd.perm"=emd.perm),
            class = "EMD_nonzero")
  
}
EMD_sigall <- function(emdall,Hur_gene,nonHur_gene, Hur, nonHur) {
  
  structure(list("emdall"=emdall, "Hur_gene"=Hur_gene,"nonHur_gene"=nonHur_gene, "Hur"=Hur, "nonHur"=nonHur),
            class = "EMD_sigall")
  
}

calculate_single<- function(data=data, condition=condition,Hur_gene=Hur_gene, binSize=0.2, nperm){
  
 require(arm)
  require(aod)
  require(fdrtool)
  
  names(condition) <- colnames(data)
  
  data<-as.matrix(data)
  databinary<- data
  databinary[databinary>0] <- 1
  
  if(is.null(Hur_gene)){
    beta_glm <- bayeswald(databinary,condition)
    
    Hur_gene <- rownames(beta_glm)[beta_glm[,"chi-test"]<0.05]
  }
 
  nonHur_gene <- setdiff(rownames(databinary),(Hur_gene))
  if(length(Hur_gene)!=0){
    message("do EMD for Hur genes.")
    Hur.result <- calculate_emd_hur(data[Hur_gene, ], condition, binSize, nperm)
  }
  
  if(length(nonHur_gene)!=0){
    message("do EMD for nonHur genes.")
    nonHur.result <- calculate_nonzero(data[nonHur_gene,], databinary[nonHur_gene,], condition, nperm, binSize)
  }
  
  if(length(Hur_gene)==0){
    emdall<- nonHur.result$emd
    Hur.result<-list()
    
  } else if (length(nonHur_gene)==0){
    emdall<- Hur.result$emd
    nonHur.result<-list()
    
  } else if (length(Hur_gene)!=0 & length(nonHur_gene)!=0){
    emdall<- rbind(Hur.result$emd, nonHur.result$emd)
    
  }
  
  emdall<- emdall[rownames(data),]
  #fdr <- fdrtool(emdall[,2], statistic="pvalue")
  padjust<-p.adjust(emdall[,2],method = "fdr", n= nrow(emdall))
  #emdall<- (cbind(emdall,padjust,fdr$qval))
  emdall<- (cbind(emdall,padjust))
  colnames(emdall)<- c("emd","pvalue","padjust")

  
  EMD_sigall(emdall,Hur_gene,nonHur_gene, Hur.result, nonHur.result)
  
}

bayeswald<- function(databinary, condition){
  cdr<-scale(colSums((databinary)>0))
  names(cdr)<- condition
  
  beta_glm <- apply(databinary[1:dim(databinary)[1],], 1, .bayeslogit, condition, cdr)
  rownames(beta_glm)<- c("intercept","condition","cdr","chi-test")
  t(beta_glm)
}

.bayeslogit<- function(x,condition, cdr){
  M3<- bayesglm (x ~ factor(condition) + cdr, family=binomial(link="logit"))
  waldt<-wald.test(b=coef(object=M3), Sigma=vcov(object=M3), Terms=2)
  c(M3$coefficients, waldt$result$chi2[3])
}

calculate_nonzero <- function(data, databinary, condition, nperm, binSize){
  ## prepare data 
 
  data.df<-list()
  outcomes<-list()
  sample_names<-list()
  datalist<-list()
  
  for (i in 1:dim(data)[1]){
    data.df[[i]] <- data[i,databinary[i,]>0]
    names(data.df[[i]])<- colnames(data)[databinary[i,]>0]
    outcomes[[i]] <- condition[databinary[i,]>0]
    sample_names[[i]] <- colnames(data)[databinary[i,]>0]
    
    if(length(setdiff(unique(condition),unique(outcomes[[i]])))!=0){
      diff <- setdiff(unique(condition),unique(outcomes[[i]]))
      for (j in diff){
        temp<-names(outcomes[[i]])
        outcomes[[i]] <- c(outcomes[[i]],j)
        names(outcomes[[i]]) <- c(temp,j)
        
        temp<-names(data.df[[i]])
        data.df[[i]] <- c(data.df[[i]],0)
        names(data.df[[i]])<- c(temp,j)
        sample_names[[i]] <- c(sample_names[[i]],j)
      }
    }
    datalist[[i]]<- list(data.df=(data.df[[i]]),outcomes=(outcomes[[i]]),sample_names=(sample_names[[i]]) ) 
  }
  names(datalist)<-rownames(data)[1:dim(data)[1]]
  
  
  ## calculate emd for each gene
 message("calculate emd for each gene...", appendLF=FALSE)
  emd <- sapply(datalist, calculate_emd_gene_sig, binSize)
  message("done")

  ## calculate emd.perm # ------------------ perm --------------------
  message("calculate emd.perm...", appendLF=FALSE)
  if(nperm==1){
    emd.perm <- (sapply(datalist, .emd_pairwise_perm_sig , nperm, binSize))
  } else {
    emd.perm <- t(sapply(datalist, .emd_pairwise_perm_sig , nperm, binSize))
  }
  message("done")
  
  ## calculate pvalue # ------------------ p-values --------------------
  message("calculate q-value...", appendLF=FALSE)
  pvals <- .emd_pairwise_qval_sig(emd,emd.perm,nperm)
  message("done")
  
  EMD<- cbind(emd, pvals)
  EMD_nonzero(data.df, outcomes, EMD, emd.perm)

}

# ------------------ calculate emd for one gene ---------------------
calculate_emd_gene_sig <- function(datalist, binSize=0.2) {

 
  vec<-datalist[[1]]
  outcomes<- datalist[[2]]
  sample_names<- datalist[[3]]
  
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  
  EMD.tab <- matrix(NA, nrow=1, ncol=dim(pairs)[2])
  
  
  for (p in 1:dim(pairs)[2])
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    EMD <- .emd_gene_pairwise(vec,src.lab,sink.lab,binSize)
    EMD.tab[1,p] <- EMD
  }
  
  EMD.tab <- as.numeric(EMD.tab)
  
  max(EMD.tab,na.rm=TRUE)
}

.adddata <- function(miss_con, vec, outcomes, sample_names){
  outcomes <- c(outcomes,miss_con)
  temp<-names(vec)
  vec <- c(vec,0)
  names(vec)<- c(temp,miss_con)
  sample_names <- c(sample_names,miss_con)
} 



.emd_pairwise_perm_sig <- function(datalist, nperm, binSize) {

  vec<-datalist[[1]]
  outcomes<- datalist[[2]]
  sample_names<- datalist[[3]]
  
  
  # ------------------ calculate permuted emd scores ---------------------
  sample_count <- length(outcomes)
  # matrix to hold permuted emd values
  emd.perm <- matrix(nrow=1, ncol=nperm)
  colnames(emd.perm) <- as.character(1:nperm)
  
  for (i in 1:nperm) {
    # permute samples
    idx.perm <- sample(1:sample_count, replace=FALSE)
    sample.id <- sample_names
    outcomes.perm <- outcomes[idx.perm]
    names(outcomes.perm) <- sample.id
    
    # calculate emd for permuted samples
    perm.val<-calculate_emd_gene(vec,outcomes.perm, sample_names, binSize)
    
    emd.perm[i] <- unlist(sapply(perm.val,"[",1))
    

  }
  emd.perm
}   

# ------------------ calculate pvalue ---------------------
.emd_pairwise_qval_sig <- function(emd, emd.perm,nperm ) { 
  
   temp<- cbind(emd,emd.perm)
   pvals<- apply(temp,1, function(x){sum(x[-1]>=as.numeric(x[1]))/nperm})

  # pvals<- apply(temp,1, function(x){abs(sum((x[-1])>= abs(as.numeric(x[1])) ))/(nperm)})
  
}

