source("helper_functions.R")

#############################################################################################
# STEP 1: preparing putatinve interaction 
# We have aggregated and process putative interctions between miRNA-RNA target and TF-RNA target in 
# Sample putative interactions and gene expression data is contained in object regression.data
# Read Readme.md for more detais 

load("sampleData.rda")
# create a regression.data list object that contains all input data needed for subsequent steps
regression.data = list(RNA = RNA, miRNA = miRNA, cna = cna, methyl = methyl, 
                       tf.target.interactions = tf.target.interactions, 
                       miRNA.target.interactions = miRNA.target.interactions)
#############################################################################################
# STEP 2: Indentify miRNA-RNA interactions 
# identify frequently selected predictors from the first 100 LASSO runs
coefs = get_nonzero_coef_lasso(regression.data = regression.data)
# get frequency of selected predictors 
bt.interval.dt = get_bootstrap_confidence_interval(regression.data = regression.data)
# find regulator-target interactions 
regulator.target.pair.dt = get_regulator_pair_list(coefs = coefs, 
                                                   bt.interval.dt = bt.interval.dt, 
                                                   RNA.targets = rownames(regression.data$RNA))
# select miRNA-RNA target interactions 
target.miRNA.list = get_target_miRNA_list(regulator.target.pair.dt)

#############################################################################################
# STEP 3: Indentify ceRNA interactions 

# create candidate ceRNA pairs (RNAs share at least one common miRNA regulator) 
pair.dt = create_pair_dt(RNA.miRNA.list = target.miRNA.list)
# compute hypergeometric p-value to test the significance of number of shared miRNAs for each pair
pair.dt = get_hypergeometric_pvalue(pair.dt = pair.dt,
                                    RNA.miRNA.list = target.miRNA.list)
# compute correlation for each pair 
pair.dt = get_correlation(pair.dt = pair.dt, RNA.expression = regression.data$RNA)
# compute partial correlation and sensitivity for each pair based on their shared RNA regulators 
pair.dt = get_partial_correlation(pair.dt = pair.dt,
                                  RNA.miRNA.list = target.miRNA.list,
                                  RNA.regression = regression.data$RNA,
                                  miRNA.regression = regression.data$miRNA)
# compute empirical pvalue of sensitivity correlation for each pair
permutated.sensitivity.cor.list = get_permutated_sensitivity_correlation_list(pair.dt = pair.dt, 
                                                                              RNA.miRNA.list = target.miRNA.list,
                                                                              RNA = regression.data$RNA,
                                                                              miRNA = regression.data$miRNA)

pair.dt = get_empirical_pvalue_sensitivity_correlation(pair.dt = pair.dt, 
                                                       permutated.sensitivity.cor.list = permutated.sensitivity.cor.list)


