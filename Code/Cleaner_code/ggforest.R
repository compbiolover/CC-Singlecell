#ggforest.R

active_predictors <-paste(my_best_genes, collapse = "+")
regular_cox_df <- cox_df[,my_best_genes]
regular_cox_df$days.to.last.follow.up <- cox_df$days.to.last.follow.up
regular_cox_df$vital.status <- cox_df$vital.status
my_formula <- paste("~", paste(active_genes[1:length(my_best_genes)], collapse = "+"))
regular_cox2 <- coxph(Surv(time = regular_cox_df$days.to.last.follow.up, event = regular_cox_df$vital.status)~., data = regular_cox_df)
