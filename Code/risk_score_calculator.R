#Name: risk_score_calculator.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Calculate Kaplan-Meier risk estimators and plot the results-----
risk_score_calculator <- function(my.file="Data/Data-from-Cleaner-code/Regular_cox_model_outputs/coad_and_read_regular_cox_genes_1800.csv", 
                                  tumor.data=FALSE, n.data=FALSE,
                                  my.title="Finished KM Plot",
                                  my.x = "Time (days)",
                                  my.y = "Survival probability",
                                  my.digits = 4,
                                  cox.df=cox_df,
                                  show.pval=TRUE,
                                  show.confidence = FALSE,
                                  show.pval.method=FALSE,
                                  my.line.size = 3.0,
                                  my.censor.size = 4){
  #Required packages----
  require(survival)
  
  #Required functions----
  risk_converter <-function(my.name=my_genes[1], my.data=risk_df, my.med.exp=gene_info){
    if(my.med.exp[counter]>0){
      my.data[,my.name]<- ifelse(my.data[,my.name]> my.med.exp[,my.name], 1, 0)
    }else{
      my.data[,my.name]<- ifelse(my.data[,my.name]< my.med.exp[,my.name], -1, 0)
    }
    return(my.data[,my.name])
  }
  km_plotter <- function(km.fit        =km_fit,
                         data.source   =surv_gene_df,
                         p.value       =show.pval,
                         pval.digits   =my.digits,
                         confidence.int=show.confidence,
                         legend.labs   =c("High risk", "Low risk"),
                         legend.title  =my.y,
                         x.lab         =my.x,
                         plot.title    =my.title,
                         color.pal     =c("red", "blue"),
                         size          = my.line.size,
                         censor.size   = my.censor.size,
                         my.km.plot    = "KM_plot.svg",
                         my.km.plot.type = "svg",
                         my.plot.height  = 34,
                         my.plot.width   = 34,
                         my.plot.dpi     = 300,
                         my.plot.x.size  = 40,
                         my.plot.y.size  = 40,
                         my.plot.legend.size = 40,
                         my.plot.main.size = 40){
    
    #Loading the needed package. If not installed it is automatically done.----
    require(survminer)
    
    custom_theme <- function() {
      theme_survminer() %+replace%
        theme(
          plot.title=element_text(hjust=0.5),
          axis.line = element_line(size = 1.5)
        )
    }
    
    #KM Curves plotting code----
    sur_Plot<-ggsurvplot(km.fit,
                         data=data.source,
                         pval=show.pval,
                         pval.method = show.pval.method,
                         pval.size=pval.digits,
                         conf.int = confidence.int,
                         legend.labs=legend.labs,
                         legend.title=legend.title,
                         xlab=x.lab,
                         title=plot.title,
                         palette=color.pal,
                         font.main= c(40, "bold"),
                         font.x= c(40, "bold"),
                         font.y=c(40, "bold"),
                         font.tickslab=c(40, "plain"),
                         font.legend=c(40, "plain"),
                         size = size,
                         censor.size=3,
                         censor.shape="|",
                         ggtheme = custom_theme())
    
    #Saving the plot to .svg or other specified format
    ggsave(filename = my.km.plot,
           plot = print(sur_Plot$plot, newpage = FALSE),
           device= my.km.plot.type, dpi=my.plot.dpi,
           width = my.plot.width, height = my.plot.height, units = "cm")
    
    #Returning our finished KM plot----
    return(sur_Plot)
  }
  
  
  #List to store the return objects in
  survival_return <- vector(mode = "list", length = 2)
  
  #Read in the data----
  my_file <- read.csv(my.file)
  my_file <- my_file[,2:3]
  colnames(my_file)[c(1,2)] <- c("gene","coef")
  
  if(tumor.data==FALSE & n.data==FALSE){
    risk_df <- cox.df[,my_file$gene]
    risk_df <- as.matrix(risk_df)
    print(dim(risk_df))
    gene_sign <- ifelse(my_file$coef>0, 1,-1)
    print(length(gene_sign))
    risk_df <- risk_df%*%diag(gene_sign)
    risk_df <- as.data.frame(risk_df)
    colnames(risk_df) <- my_file$gene
    gene_info <- data.frame(med_expression=apply(risk_df, 2, median))
    gene_info <- t(gene_info)
    gene_info <- as.data.frame(gene_info)
    my_genes <- colnames(risk_df[1:length(colnames(risk_df))])
    risk_df$vital.status <- cox.df$vital.status
    risk_df$time <- cox.df$days.to.last.follow.up
    counter <- 1
    my_converted_scores <- list()
    total_risk_length <- length(colnames(risk_df))
    risk_length <- total_risk_length-2
    for(x in my_genes){
      current_risk <- risk_converter(my.name = my_genes[counter],
                                     my.data = risk_df[1:risk_length],
                                     my.med.exp = gene_info[1:risk_length])
      my_converted_scores[[as.character(my_genes[counter])]] <- current_risk
      counter <- counter + 1
    }
    #Making a data frame of the converted scores
    converted_df <- data.frame(my_df=1:dim(cox.df)[1])
    counter <- 1
    for (x in my_converted_scores) {
      converted_df <- cbind(my_converted_scores[my_genes[counter]], converted_df)
      counter <- counter + 1
    }
    
    converted_df[,"my_df"] <- NULL
    rownames(converted_df) <- rownames(cox.df)
    converted_df$vital.status <- cox.df$vital.status
    converted_df$time <- cox.df$days.to.last.follow.up
    converted_df <- apply(converted_df, 2, as.numeric)
    patient_risks <- rowSums(x=converted_df[,1:risk_length])
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- patient_risks
    converted_df <- apply(converted_df, 2, as.numeric)
    median_risk <- median(abs(converted_df[,1:risk_length]))
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- ifelse(converted_df$risk>median_risk, "high", "low")
    write.csv(converted_df, "~/Desktop/km_plotter_df.csv")
    km_fit <- survfit(Surv(time, vital.status) ~ risk, data = converted_df)
    saveRDS(km_fit, "~/Desktop/km_fit_km_plotter.rds")
    finished_plot <- km_plotter(km.fit = km_fit, data.source = converted_df,
                                p.value = TRUE,
                                plot.title = my.title)
    finished_plot
    
    survival_return[["KM Plot"]] <- finished_plot
  }
  
  #Only for if there is tumor data included----
  if(tumor.data==TRUE & n.data==FALSE){
    tumor_info <- filter(my_file, gene=="tumor.stage1" | gene=="tumor.stage2" | gene=="tumor.stage3" | gene=="tumor.stage4")
    tumor_names <- tumor_info$gene
    print(tumor_names)
    tumor_info <- length(tumor_info$gene)
    just_genes <- length(my_file$gene) - tumor_info
    if(just_genes==1){
      risk_df <- cox.df[,my_file$gene[1]]
      risk_df <- as.data.frame(risk_df)
    }else{
      risk_df <- cox.df[,my_file$gene[1:just_genes]]
    }
    
    
    for(x in tumor_names){
      if(x=="tumor.stage1"){
        risk_df$tumorstage1 <- ifelse(cox.df$tumor.stage==1, 5,0)
      }else if(x=="tumor.stage2"){
        risk_df$tumorstage2 <- ifelse(cox.df$tumor.stage==2, 10,0)
      }else if(x=="tumor.stage3"){
        risk_df$tumorstage3 <- ifelse(cox.df$tumor.stage==3, 15,0)
      }else if(x=="tumor.stage4"){
        risk_df$tumorstage4 <- ifelse(cox.df$tumor.stage==4, 30,0)
      }
    }
    risk_df <- as.matrix(risk_df)
    gene_sign <- ifelse(my_file$coef>0, 1, -1)
    risk_df <- risk_df%*%diag(gene_sign)
    risk_df <- as.data.frame(risk_df)
    colnames(risk_df) <- my_file$gene
    gene_info <- data.frame(med_expression=apply(risk_df, 2, median))
    gene_info <- t(gene_info)
    gene_info <- as.data.frame(gene_info)
    my_genes <- colnames(risk_df[1:length(colnames(risk_df))])
    risk_df$vital.status <- cox.df$vital.status
    risk_df$time <- cox.df$days.to.last.follow.up
    counter <- 1
    my_converted_scores <- list()
    total_risk_length <- length(colnames(risk_df))
    risk_length <- total_risk_length-2
    for(x in my_genes){
      current_risk <- risk_converter(my.name = my_genes[counter], my.data = risk_df[1:risk_length], my.med.exp = gene_info[1:risk_length])
      my_converted_scores[[as.character(my_genes[counter])]] <- current_risk
      counter <- counter + 1
    }
    #Making a data frame of the converted scores
    converted_df <- data.frame(my_df=1:dim(cox.df)[1])
    counter <- 1
    for (x in my_converted_scores) {
      converted_df <- cbind(my_converted_scores[my_genes[counter]], converted_df)
      counter <- counter + 1
    }
    
    converted_df[,"my_df"] <- NULL
    rownames(converted_df) <- rownames(cox.df)
    converted_df$vital.status <- cox.df$vital.status
    converted_df$time <- cox.df$days.to.last.follow.up
    converted_df <- apply(converted_df, 2, as.numeric)
    patient_risks <- rowSums(x=converted_df[,1:risk_length])
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- patient_risks
    converted_df <- apply(converted_df, 2, as.numeric)
    median_risk <- median(abs(converted_df[,1:risk_length]))
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- ifelse(converted_df$risk>median_risk, "high", "low")
    km_fit <- survfit(Surv(time, vital.status) ~ risk, data = converted_df)
    finished_plot <-km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = my.title)
    finished_plot
    survival_return[["KM Plot"]] <- finished_plot
    
  }else if(tumor.data==TRUE & n.data==TRUE){
    tumor_info <- filter(my_file, gene=="tumor.stage1" | gene=="tumor.stage2" | gene=="tumor.stage3" | gene=="tumor.stage4")
    n_info <- filter(my_file, gene=="ajcc.n1" | gene=="ajcc.n2" | gene=="ajcc.n3")
    tumor_names <- tumor_info$gene
    n_names <- n_info$gene
    print(tumor_names)
    print(n_names)
    tumor_info <- length(tumor_info$gene)
    n_info <- length(n_info$gene)
    just_genes <- length(my_file$gene) - tumor_info - n_info
    if(just_genes==1){
      risk_df <- cox.df[,my_file$gene[1]]
      risk_df <- as.data.frame(risk_df)
    }else{
      risk_df <- cox.df[,my_file$gene[1:just_genes]]
    }
    
    for(x in tumor_names){
      if(x=="tumor.stage1"){
        risk_df$tumorstage1 <- ifelse(cox.df$tumor.stage==1, 5,0)
      }else if(x=="tumor.stage2"){
        risk_df$tumorstage2 <- ifelse(cox.df$tumor.stage==2, 10,0)
      }else if(x=="tumor.stage3"){
        risk_df$tumorstage3 <- ifelse(cox.df$tumor.stage==3, 15,0)
      }else if(x=="tumor.stage4"){
        risk_df$tumorstage4 <- ifelse(cox.df$tumor.stage==4, 30,0)
      }
    }
    
    for(x in n_names){
      if(x=="ajcc.n0"){
        risk_df$ajcc.n0 <- ifelse(cox.df$ajcc.n==0, 5,0)
      }else if(x=="ajcc.n1"){
        risk_df$ajcc.n1 <- ifelse(cox.df$ajcc.n==1, 10,0)
      }else if(x=="ajcc.n2"){
        risk_df$ajcc.n2 <- ifelse(cox.df$ajcc.n==2, 15,0)
      }else if(x=="ajcc.n3"){
        risk_df$ajcc.n3 <- ifelse(cox.df$ajcc.n==3, 20,0)
      }
    }
    
    
    risk_df <- as.matrix(risk_df)
    gene_sign <- ifelse(my_file$coef>0, 1, -1)
    print(length(gene_sign))
    print(dim(risk_df))
    risk_df <- risk_df%*%diag(gene_sign)
    risk_df <- as.data.frame(risk_df)
    colnames(risk_df) <- my_file$gene
    gene_info <- data.frame(med_expression=apply(risk_df, 2, median))
    gene_info <- t(gene_info)
    gene_info <- as.data.frame(gene_info)
    my_genes <- colnames(risk_df[1:length(colnames(risk_df))])
    risk_df$vital.status <- cox.df$vital.status
    risk_df$time <- cox.df$days.to.last.follow.up
    counter <- 1
    my_converted_scores <- list()
    total_risk_length <- length(colnames(risk_df))
    risk_length <- total_risk_length-2
    for(x in my_genes){
      current_risk <- risk_converter(my.name = my_genes[counter], my.data = risk_df[1:risk_length], my.med.exp = gene_info[1:risk_length])
      my_converted_scores[[as.character(my_genes[counter])]] <- current_risk
      counter <- counter + 1
    }
    #Making a data frame of the converted scores
    converted_df <- data.frame(my_df=1:dim(cox.df)[1])
    counter <- 1
    for (x in my_converted_scores) {
      converted_df <- cbind(my_converted_scores[my_genes[counter]], converted_df)
      counter <- counter + 1
    }
    
    converted_df[,"my_df"] <- NULL
    rownames(converted_df) <- rownames(cox.df)
    converted_df$vital.status <- cox.df$vital.status
    converted_df$time <- cox.df$days.to.last.follow.up
    converted_df <- apply(converted_df, 2, as.numeric)
    patient_risks <- rowSums(x=converted_df[,1:risk_length])
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- patient_risks
    converted_df <- apply(converted_df, 2, as.numeric)
    median_risk <- median(abs(converted_df[,1:risk_length]))
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- ifelse(converted_df$risk>median_risk, "high", "low")
    km_fit <- survfit(Surv(time, vital.status) ~ risk, data = converted_df)
    finished_plot <-km_plotter(km.fit = km_fit, data.source = converted_df, p.value = TRUE, plot.title = my.title)
    survival_return[["KM Plot"]] <- finished_plot
  }
  
  survival_return[["Survival DF"]] <- converted_df
  return(survival_return)
}


