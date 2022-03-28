#Name: miRNA_heatmap.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Allows me to do background running of this code

#The loop for miRNA and miRNA target generation----
mirna_targets <- seq(10, 1010, by=100)
mirna_num <- 900

for(mt in mirna_targets[1:11]){
  mirna.genes <- mirna_calculator(cancer.type1 = "colon cancer",
                                  cancer.type2 = "rectal cancer",
                                  max.miR.targets = mt,
                                  cancer.up = TRUE,
                                  mirna.remove = c("hsa-miR-129-2-3p",
                                                   "hsa-miR-129-1-3p",
                                                   "hsa-miR-454-3p",
                                                   "hsa-miR-365a-5p",
                                                   "hsa-miR-873-5p",
                                                   "hsa-miR-505-3p",
                                                   "hsa-miR-496",
                                                   "hsa-miR-483-3p",
                                                   "hsa-miR-455-3p",
                                                   "hsa-miR-411-5p",
                                                   "hsa-miR-203a-3p",
                                                   "hsa-miR-183-5p",
                                                   "hsa-miR-133a-3p"
                                                   ),
                                  max.mirnas = mirna_num,
                                  save.mirna.genes = TRUE,
                                  mirna.ranking.name = paste0("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/Data/TCGA-READ/Mirna/Global_Search/Inputs/RData/",mirna_num,"_",mt,"_targets_mirna_ranking.rds"))

}


#Just for testing whether dbDEMC only miRNAs make a difference compared to the
#overlap of miRMap and dbDEMC
# mirna_targets <- seq(10,1010,100)
# for(mt in mirna_targets){
#   mirna.genes <- mirna_calculator_dbdemc_only(cancer.type1 = "colon cancer",
#                                               cancer.type2 = "rectal cancer",
#                                               max.miR.targets = mt,
#                                               cancer.up = TRUE,
#                                               mirna.remove = c("hsa-miR-129-2-3p",
#                                                                "hsa-miR-129-1-3p",
#                                                                "hsa-miR-454-3p",
#                                                                "hsa-miR-365a-5p",
#                                                                "hsa-miR-873-5p",
#                                                                "hsa-miR-505-3p",
#                                                                "hsa-miR-496",
#                                                                "hsa-miR-483-3p",
#                                                                "hsa-miR-455-3p",
#                                                                "hsa-miR-411-5p",
#                                                                "hsa-miR-203a-3p",
#                                                                "hsa-miR-183-5p",
#                                                                "hsa-miR-133a-3p",
#                                                                "hsa-miR-142-3p",
#                                                                "hsa-miR-504-5p",
#                                                                "hsa-miR-126-3p"
#                                               ),
#                                               max.mirnas = 5,
#                                               save.mirna.genes = TRUE,
#                                               mirna.ranking.name = paste0("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/Data/TCGA-COAD/Mirna/Global_Search/Inputs/RData/",300,"_",mt,"_targets_dbdemc_mirna_ranking_test.rds"))
#   
#   
#   
#   
# }
# mirna.genes <- mirna_calculator_dbdemc_only(cancer.type1 = "colon cancer",
#                                             cancer.type2 = "rectal cancer",
#                                             max.miR.targets = 10,
#                                             cancer.up = TRUE,
#                                             mirna.remove = c("hsa-miR-129-2-3p",
#                                                              "hsa-miR-129-1-3p",
#                                                              "hsa-miR-454-3p",
#                                                              "hsa-miR-365a-5p",
#                                                              "hsa-miR-873-5p",
#                                                              "hsa-miR-505-3p",
#                                                              "hsa-miR-496",
#                                                              "hsa-miR-483-3p",
#                                                              "hsa-miR-455-3p",
#                                                              "hsa-miR-411-5p",
#                                                              "hsa-miR-203a-3p",
#                                                              "hsa-miR-183-5p",
#                                                              "hsa-miR-133a-3p"
#                                                              ),
#                                             max.mirnas = 100,
#                                             save.mirna.genes = TRUE,
#                                             mirna.ranking.name = paste0("~/Documents/PhD Program/Hong Lab/Projects/CC_Singlecell/Data/TCGA-READ/Mirna/Global_Search/Inputs/RData/",100,"_",10,"_targets_dbdemc_mirna_ranking.rds"))
