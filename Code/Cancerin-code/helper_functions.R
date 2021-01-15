library(data.table)
library(pbapply)
library(doParallel)
library(parallel)
library(HDCI)
library(glmnet)
library(bnlearn)

# get regression df for each RNA with cna, methyl, TF and miRNA included as candidate regulators
getDF = function(regression.data, RNA, include.lncRNA = F){
  miRNA.target.interactions = regression.data$miRNA.target.interactions
  miRNA.candidates = as.character(miRNA.target.interactions[miRNA.target.interactions$target == RNA,c("miRNA")])
  miRNA.candidates = miRNA.candidates[miRNA.candidates %in% rownames(regression.data$miRNA)]
  
  tf.target.interactions = regression.data$tf.target.interactions
  TF.candidates = tf.target.interactions$tf[which(tf.target.interactions$target == RNA)] 
  TF.candidates = setdiff(TF.candidates[TF.candidates %in% rownames(regression.data$RNA)], RNA)
  
  rna.expression.vector = regression.data$RNA[RNA,]
  rna.cna.vector = regression.data$cna[RNA,]
  rna.methyl.vector = regression.data$methyl[RNA,]
  predicted.TF.df = regression.data$RNA[TF.candidates,]
  predicted.miRNA.df = regression.data$miRNA[miRNA.candidates,]
  
  
  if (include.lncRNA){
    if (nrow(predicted.miRNA.df) < 1){
      return(NULL)
    }else{
      regression.df = as.data.frame(scale(cbind(t(rna.expression.vector), t(predicted.miRNA.df))))
      colnames(regression.df)[1] = c("RNA")
      return(regression.df)
    }
  }
  
  if (nrow(predicted.TF.df) >= 1 & nrow(predicted.miRNA.df) < 1){
    regression.df = as.data.frame(scale(cbind(t(rna.expression.vector),t(rna.cna.vector),t(rna.methyl.vector),
                                              t(predicted.TF.df))))
  }else if (nrow(predicted.TF.df) < 1 & nrow(predicted.miRNA.df) >= 1){
    regression.df = as.data.frame(scale(cbind(t(rna.expression.vector),t(rna.cna.vector),t(rna.methyl.vector),
                                              t(predicted.miRNA.df))))
  }else if (nrow(predicted.TF.df) >= 1 & nrow(predicted.miRNA.df) >= 1){
    regression.df = as.data.frame(scale(cbind(t(rna.expression.vector),t(rna.cna.vector),t(rna.methyl.vector),
                                              t(predicted.miRNA.df),t(predicted.TF.df))))
  }else{
    regression.df = NULL
  }
  
  if (is.null(regression.df)){
    return(NULL)
  }
  
  colnames(regression.df)[1:3] = c("RNA","CNA","Methyl")
  regression.df$CNA[which(is.na(regression.df$CNA))] = 0
  regression.df$Methyl[which(is.na(regression.df$Methyl))] = 0
  
  return(regression.df)
}

get_nonzero_coef_lasso = function(regression.data){
  RNA.targets = rownames(regression.data$RNA)
  outfile = paste("log/", cancer.type ,".normal.out",sep = "")
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster, 
                          varlist = c("RNA.targets","regression.data", "getDF",
                                      "perform_boot_Lasso"),
                          envir = environment())
  invisible(clusterEvalQ(cl = cluster, library(data.table))) 
  invisible(clusterEvalQ(cl = cluster, library(glmnet))) 
  invisible(clusterEvalQ(cl = cluster, library(HDCI)))
  
  cat("Start variable selection: normal lasso ... \n")
  
  coefs = pbapply::pblapply(cl = cluster, X = 1:length(RNA.targets), FUN = function(RNA.index){
    RNA = RNA.targets[RNA.index]
    cat("Start", RNA, "\n")
    include.lncRNA = grepl(RNA, pattern = "lncRNA")
    regression.df = getDF(regression.data = regression.data, RNA = RNA, include.lncRNA = include.lncRNA)
    if (is.null(regression.df)){
      return(NULL)
    }
    
    y.matrix = as.matrix(regression.df$RNA)
    x.matrix = as.matrix(regression.df[,-1])
    colnames(x.matrix) = colnames(regression.df)[-1]
    candidate.regulators = colnames(x.matrix)
    # if only 1 miRNA
    if(ncol(x.matrix) == 1){
      if(cor(x.matrix, y.matrix) > 0){
        return(NULL)
      }else{
        if (cor.test(x.matrix, y.matrix)$p.value >= 0.05){
          return(NULL)
        }else{
          lasso.dt = data.frame(regulator=candidate.regulators, target=RNA, 
                                coef = cor(x.matrix, y.matrix), 
                                stringsAsFactors = F)
          lasso.dt = as.data.table(lasso.dt)
          lasso.dt$pair = paste(lasso.dt$regulator, lasso.dt$target,sep="~")
          return(lasso.dt)
        }
      }
    }
    
    # get selected coefs
    selected.coefs = lapply(1:100, function(iter){
      #cat(iter, "\n")
      y.matrix = as.matrix(regression.df$RNA)
      x.matrix = as.matrix(regression.df[,-1])
      colnames(x.matrix) = colnames(regression.df)[-1]
      
      lasso.result = tryCatch({
        HDCI::Lasso(x = x.matrix, y = y.matrix,
                    fix.lambda = F,nfolds = 10,
                    cv.method = "cv1se",
                    standardize = T, parallel = F)
      }, warning = function(w) {
        NULL
      }, error = function(e) {
        NULL
      })
      
      if (is.null(lasso.result)){
        return(NULL)
      }
      
      coefs = lasso.result$beta
      names(coefs) =  colnames(x.matrix) 
      coefs = coefs[which(coefs != 0)]
      return(coefs)
    })
    
    # make consistently selected data structure
    lasso.dt = NULL
    dt.list = lapply(1:length(selected.coefs), function(target.index){
      regulators = selected.coefs[[target.index]]
      if (is.null(regulators)){
        return(NULL)
      }
      target = RNA
      regulator.names = names(regulators)
      regulator.coef =unname(regulators)
      if (length(regulator.names) > 0){
        df = data.frame(regulator=regulator.names,
                        target = rep(target,length(regulator.names)),
                        coef=as.vector(regulator.coef),
                        stringsAsFactors = F)
        df$pair = paste(df$regulator,df$target,sep="~")
        dt = as.data.table(df)
        return(dt)
      }else{
        return(NULL)
      }
    })
    if (is.null(dt.list)){
      lasso.dt = NULL
      return(lasso.dt)
    }
    
    lasso.dt = tryCatch({
      rbindlist(dt.list)
    }, warning = function(w) {
      NULL
    }, error = function(e) {
      NULL
    })
    
    cat("Done", RNA, "\n")
    return(lasso.dt)
  })
  stopCluster(cluster)
  names(coefs) = RNA.targets
  return(coefs)
}

get_bootstrap_confidence_interval = function(regression.data){
  RNA.targets = rownames(regression.data$RNA)
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster, 
                          varlist = c("RNA.targets","regression.data", "getDF",
                                      "perform_boot_Lasso"),
                          envir = environment())
  invisible(clusterEvalQ(cl = cluster, library(data.table))) 
  invisible(clusterEvalQ(cl = cluster, library(glmnet))) 
  invisible(clusterEvalQ(cl = cluster, library(HDCI)))
  
  bt.interval.list = pbapply::pblapply(cl = cluster, 1:length(RNA.targets),FUN = function(RNA.index){
    RNA = RNA.targets[RNA.index]
    include.lncRNA = grepl(RNA, pattern = "lncRNA")
    regression.df = getDF(regression.data = regression.data, RNA = RNA, include.lncRNA = include.lncRNA)
    
    if (is.null(regression.df)){
      return(NULL)
    }
    
    y.matrix = as.matrix(regression.df$RNA)
    x.matrix = as.matrix(regression.df[,-1])
    colnames(x.matrix) = colnames(regression.df)[-1]
    
    if(ncol(x.matrix) == 1){
      if(cor(x.matrix, y.matrix) > 0){
        return(NULL)
      }else{
        if (cor.test(x.matrix, y.matrix)$p.value >= 0.05){
          return(NULL)
        }else{
          bt.interval = matrix(c(-1000,1000), ncol = 1)
          colnames(bt.interval) = colnames(x.matrix)
          return(bt.interval)
        }
      }
    }
    
    num.bt.replications = 500
    btlasso.result = tryCatch({
      perform_boot_Lasso(x = x.matrix, y = y.matrix,
                         B = num.bt.replications,
                         cv.method = "cv1se",
                         standardize = T)
    }, warning = function(w) {
      NULL
    }, error = function(e) {
      NULL
    })
    
    if (is.null(btlasso.result)){
      return(NULL)
    }
    
    bt.interval = btlasso.result$interval
    colnames(bt.interval) = colnames(x.matrix)
    return(bt.interval)
  })
  stopCluster(cluster); gc()
  names(bt.interval.list) = RNA.targets
  return(bt.interval.list)
}

# bootstrap LASSO 
perform_boot_Lasso = function (x, y, B = 500, type.boot = "residual", alpha = 0.05, 
                               cv.method = "cv", nfolds = 10, foldid, cv.OLS = FALSE, tau = 0, 
                               standardize = TRUE, intercept = TRUE){
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  selectset <- rep(0, p)
  Beta <- rep(0, p)
  globalfit <- glmnet(x, y, standardize = standardize, intercept = intercept)
  lambda <- globalfit$lambda
  cvfit <- escv.glmnet(x, y, lambda = lambda, nfolds = nfolds, 
                       tau = tau, cv.OLS = cv.OLS, parallel = F, standardize = standardize, 
                       intercept = intercept)
  if (cv.method == "cv") {
    lambda.opt <- cvfit$lambda.cv
  }
  if (cv.method == "escv") {
    lambda.opt <- cvfit$lambda.escv
  }
  if (cv.method == "cv1se") {
    lambda.opt <- cvfit$lambda.cv1se
  }
  fitlasso <- predict(globalfit, type = "coefficients", s = lambda.opt)
  fit.value <- predict(globalfit, newx = x, s = lambda.opt)
  betalasso <- fitlasso[-1]
  Beta <- betalasso
  
  Beta.boot <- matrix(0, nrow = B, ncol = p)
  out <- list()
  
  for (i in 1:B) {
    resam <- sample(1:n, n, replace = TRUE)
    rx <- x[resam, ]
    ry <- y[resam]
    boot.obj <- Lasso(x = rx, y = ry, lambda = lambda.opt, 
                      standardize = standardize, intercept = intercept)
    out[[i]] <- boot.obj$beta
  }
  
  
  for (i in 1:B) {
    Beta.boot[i, ] <- out[[i]]
    out[[i]] <- 0
  }
  out <- NULL
  interval <- ci(Beta, Beta.boot, alpha = alpha, type = "quantile")
  
  object <- list(lambda.opt = lambda.opt, Beta = Beta, interval = interval)
  object
}

get_regulator_pair_list = function(coefs, bt.interval.dt, RNA.targets){
  regulator.list = lapply(1:length(coefs), function(index){
    gene.name = RNA.targets[index]
    coef.dt = coefs[[gene.name]]
    bt.interval.dt = bt.interval.dt[[gene.name]]
    list(coef.dt=coef.dt, bt.interval.dt=bt.interval.dt)
  })
  names(regulator.list) = RNA.targets
  
  regulator.target.pair.list = lapply(1:length(regulator.list), function(gene.index) {
    gene.pair.info.list = regulator.list[[gene.index]]
    gene.coef.dt = gene.pair.info.list$coef.dt
    if (nrow(gene.coef.dt) == 0 || is.null(gene.coef.dt)) {
      return(NULL)
    }
    gene.coef.stat.dt = gene.coef.dt[, list(count = .N, 
                                            median.coef = median(.SD$coef)),
                                     by = "pair"]
    gene.coef.stat.dt$regulator = sapply(strsplit(gene.coef.stat.dt$pair, split = "[~]"), function(v)
      v[1])
    gene.coef.stat.dt$target = sapply(strsplit(gene.coef.stat.dt$pair, split = "[~]"), function(v)
      v[2])
    gene.coef.stat.dt = gene.coef.stat.dt[, c(4, 5, 1, 2, 3), ]
    gene.bt.interval.dt = tryCatch({
      gene.pair.info.list$bt.interval.dt
    }, warning = function(w) {
      NULL
    }, error = function(e) {
      NULL
    })
    
    if (is.null(gene.bt.interval.dt) || nrow(gene.bt.interval.dt) == 0){
      return(NULL)
    }
    
    lower.percentile = sapply(1:nrow(gene.coef.stat.dt), function(row.index) {
      regulator = gene.coef.stat.dt$regulator[row.index]
      return(gene.bt.interval.dt[1, regulator])
    })
    gene.coef.stat.dt$lower.percentile = lower.percentile
    
    upper.percentile = sapply(1:nrow(gene.coef.stat.dt), function(row.index) {
      regulator = gene.coef.stat.dt$regulator[row.index]
      return(gene.bt.interval.dt[2, regulator])
    })
    gene.coef.stat.dt$upper.percentile = upper.percentile
    
    gene.coef.stat.dt$confidence =  (gene.coef.stat.dt$median.coef >= lower.percentile) &
      (gene.coef.stat.dt$median.coef <= upper.percentile) & 
      (lower.percentile*upper.percentile > 0)
    
    return(gene.coef.stat.dt)
  })
  regulator.target.pair.dt = rbindlist(regulator.target.pair.list)
  return(regulator.target.pair.dt)
}


get_target_miRNA_list = function(regulator.target.pair.dt){
  RNA.miRNA.dt = regulator.target.pair.dt[grep(regulator, pattern = "hsa")]
  RNA.miRNA.dt = RNA.miRNA.dt[which(RNA.miRNA.dt$confidence == T & RNA.miRNA.dt$count >=75 &  RNA.miRNA.dt$median.coef < 0)]
  target.miRNA.list = lapply(unique(RNA.miRNA.dt$target), function(target){
    unique(RNA.miRNA.dt$regulator[which(RNA.miRNA.dt$target == target)])
  })
  names(target.miRNA.list) = unique(RNA.miRNA.dt$target)
  target.miRNA.list
}

# create candidate ceRNA pair from list
create_pair_dt = function(RNA.miRNA.list){
  pair.dt = as.data.table(t(combn(names(RNA.miRNA.list),2)))
  colnames(pair.dt) = c("RNAi", "RNAj"); setkeyv(pair.dt, c("RNAi","RNAj"))
  pair.dt$num.common.miRNAs = apply(pair.dt, 1, function(pair){
    RNAi = pair[1]
    RNAj = pair[2]
    num.common.miRs = length(intersect(RNA.miRNA.list[[RNAi]],RNA.miRNA.list[[RNAj]]))
    return(num.common.miRs)
  })
  pair.dt = pair.dt[pair.dt$num.common.miRNAs>0]
  return(pair.dt)
}

# Input: table of RNA pairs such that the RNAs in a pair share at least one miRNA
# Output: compute correlation and correlation's p-value for each RNA pairs
get_correlation = function(pair.dt, RNA.expression){
  cor.result = WGCNA::corAndPvalue(t(RNA.expression),use = "pairwise.complete.obs")
  corMatrix = cor.result$cor
  corPvalueMatrix = cor.result$p
  
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster,varlist = c("pair.dt", "corMatrix", "corPvalueMatrix"),
                          envir = environment())
  
  pair.dt$cor = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
    return(corMatrix[RNAi,RNAj])
  }) 
  
  pair.dt$cor.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
    return(corPvalueMatrix[RNAi,RNAj])
  }); stopCluster(cluster)
  
  return(pair.dt)
}

# Input: table of RNA pairs such that the RNAs in a pair share at least one miRNA
# Output: compute partial correlation and partial correlation's p-value for each RNA pairs
get_partial_correlation = function(RNA.regression, miRNA.regression, RNA.miRNA.list, pair.dt){
  
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster, 
                          varlist = c("RNA.regression","miRNA.regression","RNA.miRNA.list", "pair.dt"),
                          envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn,stats)))
  
  
  time = proc.time()
  pc.list = pbapply::pblapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]; ceRNAs = c(RNAi, RNAj)
    common.miRs = intersect(RNA.miRNA.list[[RNAi]], RNA.miRNA.list[[RNAj]])
    
    #RNA.RNA.df = getDataForConditionalTest(regression.data = regression.data, ceRNAs = ceRNAs, common.miRs = common.miRs)
    RNA.df = RNA.regression[ceRNAs,] # since RNA.name has only one element
    miRNA.df = miRNA.regression[common.miRs,]
    RNA.RNA.df = as.data.frame(scale(cbind(t(RNA.df), t(miRNA.df))))
    
    ci.test.result = ci.test(x = colnames(RNA.RNA.df)[1], y = colnames(RNA.RNA.df)[2], z = colnames(RNA.RNA.df)[3:ncol(RNA.RNA.df)],
                             data = RNA.RNA.df,test = "cor")
    
    return(list(pcor = ci.test.result$statistic,
                pvalue = ci.test.result$p.value))
  }); stopCluster(cluster)
  
  # get partial correlation and p.value for the test
  pair.dt$partial.corr = sapply(pc.list, function(result) result$pcor)
  pair.dt$sensitivity.cor = pair.dt$cor - pair.dt$partial.corr
  pair.dt$partial.cor.pvalue = sapply(pc.list, function(result) result$pvalue);
  return(pair.dt)
}

# Input: table of RNA pairs such that the RNAs in a pair share at least one miRNA
# Output: compute hypergeometric p-value for each RNA pairs
get_hypergeometric_pvalue = function(pair.dt, RNA.miRNA.list){
  pop.size = length(unique(unlist(RNA.miRNA.list))) # total number of miRNAs
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster,
                          varlist = c("RNA.miRNA.list", "pair.dt","pop.size","getPvalueHypergeometric"),
                          envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))
  
  pair.dt$hypergeometric.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
    miRs.RNAi = RNA.miRNA.list[[RNAi]]
    miRs.RNAj = RNA.miRNA.list[[RNAj]]
    success.pop.size = length(miRs.RNAi)
    sample.size = length(miRs.RNAj)
    success.sample.size = pair.dt$num.common.miRNAs[index]
    getPvalueHypergeometric(pop.size = pop.size, success.pop.size = success.pop.size,
                            sample.size = sample.size, success.sample.size = success.sample.size)
  })
  stopCluster(cluster)
  return(pair.dt)
}

# compute pvalue for over-enrichment test using hypergeometric test
getPvalueHypergeometric = function(pop.size, success.pop.size, sample.size, success.sample.size){
  return(phyper(success.sample.size - 1, success.pop.size, pop.size - success.pop.size, sample.size, lower.tail = F))
}

# Input: table of RNA pairs such that the RNAs in a pair share at least one miRNA
# Output: compute emperical p-value for each RNA pair's sensitivity correlation
get_permutated_sensitivity_correlation_list = function(pair.dt, RNA.miRNA.list, RNA, miRNA){
  cluster = parallel::makeCluster(detectCores() - 1)
  parallel::clusterExport(cl = cluster, varlist = c("pair.dt", "RNA", "miRNA", "getDataWithResampledMiRNAs", "ci.test"))
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn)))
  
  permutated.sensitivity.cor.list = pbapply::pblapply(X = 1:nrow(pair.dt), FUN = function(row.index){
    RNAi = pair.dt$RNAi[row.index];  RNAj = pair.dt$RNAj[row.index]
    
    ceRNA.sensitivity.permuted = sapply(1:1000, function(bootstrap_iter){
      data.expression = getDataWithResampledMiRNAs(RNAi = RNAi, RNAj = RNAj,
                                                   num.common.miRNAs = pair.dt$num.common.miRNAs[row.index],
                                                   RNA =  RNA, 
                                                   miRNA = miRNA)
      # compute permutated correlation
      ceRNA.cor = as.numeric(cor(data.expression[,1], data.expression[,2]))
      ceRNA.partial.cor = ci.test(x = colnames(data.expression)[1],
                                  y = colnames(data.expression)[2],
                                  z = colnames(data.expression)[3:ncol(data.expression)],
                                  data = data.expression,test = "cor")$statistic
      return(unname(ceRNA.cor - ceRNA.partial.cor))
    })
  })
  stopCluster(cluster)
  return(permutated.sensitivity.cor.list)
}

get_empirical_pvalue_sensitivity_correlation = function(pair.dt, permutated.sensitivity.cor.list){
  cluster = parallel::makeCluster(detectCores() - 1)
  parallel::clusterExport(cl = cluster, varlist = c("permutated.sensitivity.cor.list","pair.dt"))
  pair.dt$sensitivity.cor.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    return(sum(permutated.sensitivity.cor.list[[index]] >= pair.dt$sensitivity.cor[index])/1000)
  })
  stopCluster(cluster)
  return(pair.dt)
}
                                      
# get expression data by resampling miRNAs without replacement
getDataWithResampledMiRNAs = function(RNAi, RNAj, num.common.miRNAs, RNA, miRNA){
  RNAi.expression = t(as.matrix(RNA[RNAi,])); colnames(RNAi.expression) = RNAi
  RNAj.expression = t(as.matrix(RNA[RNAj,])); colnames(RNAj.expression) = RNAj
  # construct data expression
  random.miRs = sample(rownames(miRNA), size = num.common.miRNAs, replace = F)
  permutated.miRs.expression = as.matrix(miRNA[random.miRs, ])
  if (length(random.miRs) == 1){
    rownames(permutated.miRs.expression) = random.miRs
  }
  permutated.miRs.expression = t(miRNA[random.miRs,])
  data.expression = as.data.frame(scale(cbind(RNAi.expression, RNAj.expression, permutated.miRs.expression)))
  return(data.expression)
}

