EMDomics <- function(data, outcomes, emd, emd.perm, pairwise.emd.table, pairwise.q.table) {
  
  structure(list("data"=data, "outcomes"=outcomes,
                 "emd"=emd, "emd.perm"=emd.perm, "pairwise.emd.table"=pairwise.emd.table,
                 "pairwise.q.table"=pairwise.q.table),
            class = "EMDomics")
  
}

calculate_emd_hur <- function(data, outcomes, binSize=0.2,
                          nperm=100, pairwise.p=FALSE, seq=FALSE,
                          quantile.norm=FALSE,
                          verbose=TRUE, parallel=FALSE) {
  
  bpparam <- BiocParallel::bpparam()
  
  if (parallel == FALSE)
    bpparam <- BiocParallel::SerialParam()
  
  
  if (seq)
  {
    data<-data*1E6
    data<-log2(data+1)
  }
  
  if (quantile.norm)
  {
    preprocessCore::normalize.quantiles(data)
  }
  
  # transpose and coerce to df (for bplapply)
  data.df <- as.data.frame(t(data))
  sample_names <- rownames(data.df)
  
  # ---------- pairwise emd table -----------
  
  # generate pairwise emd table for each gene
  if (verbose)
    message("Calculating pairwise emd scores...", appendLF=FALSE)
  
  # all possible pairwise comparisons
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  names <- apply(pairs,2,function(x){paste(x[1],'vs',x[2])})
  
  emd.tab <- BiocParallel::bplapply(data.df, .emd_pairwise_table, sample_names, 
                                    outcomes, pairs,
                                    binSize, verbose,
                                    BPPARAM = bpparam)
  
  # Remove genes that were not binned properly by the histogram function
  lengths <- lapply(emd.tab, function(x){length(x)})
  remove.genes <- lengths < ncol(pairs)
  bad.lengths <- lengths[remove.genes]
  for (b in names(bad.lengths)) {
    msg <- paste('Data for gene', b, 'has been removed because it could not be binned properly.')
    message(msg)
  }
  emd.tab <- emd.tab[!remove.genes]
  data.df <- data.df[,!remove.genes]
  
  emd.tab <- matrix(unlist(emd.tab), nrow=ncol(data.df), ncol=ncol(pairs), byrow=TRUE)
  rownames(emd.tab) <- colnames(data.df)
  colnames(emd.tab) <- names
  
  if (verbose)
    message("done.")
  
  # ---------- emd ------------
  
  # calculate emd for each gene
  if (verbose)
    message("Calculating emd...", appendLF=FALSE)
  
  emd <- apply(emd.tab, 1, function(x){max(as.numeric(x),na.rm=TRUE)})
  
  emd <- as.matrix(emd)
  colnames(emd) <- "emd"
  
  if (verbose)
    message("done.")
  
  # pairwise q-value computation if specified
  if (pairwise.p) {
    emd.pairwise.q <- .emd_pairwise_q(data.df, emd, outcomes, nperm, verbose, binSize, bpparam)
  } else {
    emd.pairwise.q <- NULL
  }
  
  # --------------- permuted emd scores ----------------
  
  sample_count <- length(outcomes)
  
  # matrix to hold permuted emd values
  emd.perm <- matrix(nrow=ncol(data.df), ncol=nperm)
  rownames(emd.perm) <- colnames(data.df)
  colnames(emd.perm) <- as.character(1:nperm)
  
  for (i in 1:nperm) {
    
    msg <- paste("Calculating permuted emd #", i, " of ",
                 nperm, "...", sep="")
    
    if (verbose)
      message(msg, appendLF=FALSE)
    
    # permute samples
    idx.perm <- sample(1:sample_count, replace=FALSE)
    sample.id <- names(outcomes)
    outcomes.perm <- outcomes[idx.perm]
    names(outcomes.perm) <- sample.id
    
    # calculate emd for permuted samples
    perm.val <- BiocParallel::bplapply(data.df, calculate_emd_gene,
                                       outcomes.perm, rownames(data.df),
                                       binSize,
                                       BPPARAM = bpparam)
    
    emd.perm[,i] <- as.numeric(unlist(sapply(perm.val,"[",1)))
    
    if (verbose)
      message("done.")
    
  }
  
  # ------------------ q-values --------------------
  
  if (verbose)
    message("Calculating q-values...", appendLF=FALSE)
  
  qvals <- matrix(1, nrow=dim(data)[1], ncol=1)
  rownames(qvals) <- rownames(emd)
  colnames(qvals) <- "qvals " 
  
  for (i in 1:dim(data)[1]){
    qvals[i,]<- sum(abs(emd.perm[i,])>= abs(as.numeric(emd[i,1])))/(dim(emd.perm)[2])
    
  }

  emd.qval <- as.matrix(qvals)
  colnames(emd.qval) <- "pvalue"
  emd.qval
  
  
  if (verbose)
    message("done.")
  
  emd <- cbind(emd, emd.qval)
  
  EMDomics(data, outcomes, emd, emd.perm, emd.tab, emd.pairwise.q)
 
}


calculate_emd_gene <- function(vec, outcomes, sample_names, binSize=0.2) {

  names(vec) <- sample_names
  
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  
  EMD.tab <- matrix(NA, nrow=1, ncol=dim(pairs)[2])
  colnames<-list()
  
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


.emd_pairwise_q <- function(data.df, emd, outcomes, nperm, verbose, binSize, bpparam) {

  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  names <- apply(pairs,2,function(x){paste(x[1],'vs',x[2])})

  q.tab <- matrix(NA, nrow=ncol(data.df), ncol=ncol(pairs))
  colnames(q.tab) <- names
  rownames(q.tab) <- colnames(data.df)

  for (p in 1:ncol(pairs)) {

    # ------------------ calculate permuted emd scores ---------------------
    sample_count <- length(outcomes)

    # matrix to hold permuted emd values
    emd.perm <- matrix(nrow=ncol(data.df), ncol=nperm)
    rownames(emd.perm) <- colnames(data.df)
    colnames(emd.perm) <- as.character(1:nperm)

    msg <- paste("Beginning pairwise q-value computatiion for",names[p])
    if (verbose)
      message(msg)

    for (i in 1:nperm) {

      msg <- paste("Calculating permuted emd #", i, " of ",
                   nperm, "...", sep="")

      if (verbose)
        message(msg, appendLF=FALSE)

      # permute samples
      idx.perm <- sample(1:sample_count, replace=FALSE)
      sample.id <- names(outcomes)
      outcomes.perm <- outcomes[idx.perm]
      names(outcomes.perm) <- sample.id

      # calculate emd for permuted samples
      perm.val <- BiocParallel::bplapply(data.df, calculate_emd_gene,
                                         outcomes.perm, rownames(data.df),
                                         binSize,
                                         BPPARAM = bpparam)

      emd.perm[,i] <- unlist(sapply(perm.val,"[",1))

      if (verbose)
        message("done.")
    }

    # ------------------ q-values --------------------

    if (verbose)
      message("Calculating pairwise q-values...", appendLF=FALSE)

    qvals <- matrix(1, nrow=dim(data)[1], ncol=1)
    rownames(qvals) <- rownames(emd)
    colnames(qvals) <- "qvals "

    for (i in 1:dim(data)[1]){
      qvals[i,]<- sum(abs(emd.perm[i,])>= abs(as.numeric(emd[i,1])))/(dim(emd.perm)[2])
      
    }

    emd.qval <- as.matrix(qvals)
    colnames(emd.qval) <- "pvalue"
    emd.qval


    q.tab[,p] <- emd.qval

    if (verbose)
      message("done.")
  }

  q.tab
}

# Creates a table of all the pairwise EMD scores for one gene
.emd_pairwise_table <- function(geneData, sample_names, outcomes, pairs, binSize, verbose) {
  
  names(geneData) <- sample_names
  
  EMD.tab <- matrix(NA, nrow=1, ncol=dim(pairs)[2])
  colnames<-list()
  
  for (p in 1:dim(pairs)[2])
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    EMD <- .emd_gene_pairwise(geneData,src.lab,sink.lab,binSize)
    EMD.tab[1,p] <- EMD
  }
  
  EMD.tab <- as.numeric(EMD.tab)
  
  EMD.tab
}

# computes pairwise EMD
.emd_gene_pairwise <- function(vec, idxA, idxB, binSize=0.2) {
  dataA <- vec[idxA]
  dataB <- vec[idxB]
  
  dataA <- as.numeric(dataA)
  dataB <- as.numeric(dataB)
  
  bins <- seq(floor(min(c(dataA, dataB))),
              ceiling(max(c(dataA, dataB))),
              by=binSize )
  
  histA <- hist(dataA, breaks=bins, plot=FALSE)
  histB <- hist(dataB, breaks=bins, plot=FALSE)
  
  densA <- as.matrix(histA$density)
  densB <- as.matrix(histB$density)
  
  emdist::emd2d(densA, densB)
}

