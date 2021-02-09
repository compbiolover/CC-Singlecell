#Author: Andrew Willems <awillems@vols.utk.edu>.
#Purpose: To create Kaplan-Meir (KM) curvs from my elastic net Cox models
#R version: 4.0.2.


#Loading needed packages----
library(survminer);packageVersion("survminer")

#Now getting the model to fit.----
#Note: For the current formula data I take the gene list of the top performing
#C-index cox model.
my_survival <- Surv(days.to.last.follow.up, vital.status)
my_formula <- paste("Surv(days.to.last.follow.up, vital.status)~",paste(current_formula_data[1:length(current_formula_data)], collapse = "+"))
my_formula <- as.formula(my_formula)
my_test <- coxph(my_formula, data = df_for_train_test_split)
km_fit <- survfit(my_formula, data = df_for_train_test_split)
km_curve <- ggsurvplot(km_fit, data = df_for_train_test_split)
