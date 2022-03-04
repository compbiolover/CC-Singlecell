#Name: km_plotter.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: To plot survival data from our model fits efficiently
km_plotter <- function(km.fit                  =km_fit,
                       data.source             =surv_gene_df,
                       p.value                 =TRUE,
                       pval.size               =5,
                       show.pval.method        = FALSE,
                       confidence.int          =FALSE,
                       legend.labs             =c("High risk", "Low risk"),
                       legend.title            ="Risk",
                       x.lab                   ="Time (days)",
                       y.lab                   ="Survival probability",
                       plot.title              ="KM Plot",
                       my.plot.x.size          = 40,
                       my.plot.x.face          = "bold",
                       my.plot.y.size          = 40,
                       my.plot.y.face          = "bold",
                       my.plot.legend.size     = 40,
                       my.plot.legend.face     = "bold",
                       my.plot.main.size       = 40,
                       my.plot.main.face       = "bold",
                       my.plot.ticks.size      = 40,
                       my.plot.ticks.size.face = "plain",
                       color.pal               =c("red", "blue"),
                       surv.curv.size          = 4.5,
                       my.censor               = TRUE,
                       my.censor.size          = 4.5,
                       my.censor.shape         = "|",
                       my.km.plot              = "KM_plot.svg",
                       my.km.plot.type         = "svg",
                       my.plot.height          = 34,
                       my.plot.width           = 34,
                       my.plot.dpi             = 300,
                       my.plot.unit            = "cm"){
  
  #Loading the needed package
  require(survminer)
  
  #Custom function to modify specific aspects of 'survminer's' default theme
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(
        plot.title=element_text(hjust=0.5),
        axis.line = element_line(size = 1.5)
      )
  }
  
  #Kaplan-Meier (KM) curve plotting
  sur_Plot<-ggsurvplot(fit = km.fit,
                       data=data.source,
                       palette=color.pal,
                       conf.int = confidence.int,
                       pval=p.value,
                       pval.method = show.pval.method,
                       pval.size=pval.size,
                       legend.labs=legend.labs,
                       legend.title=legend.title,
                       xlab=x.lab,
                       ylab=y.lab,
                       title=plot.title,
                       font.main= c(my.plot.main.size, my.plot.main.face),
                       font.x= c(my.plot.x.size, my.plot.x.face),
                       font.y=c(my.plot.y.size, my.plot.y.face),
                       font.tickslab=c(my.plot.ticks.size,
                                       my.plot.ticks.size.face),
                       font.legend=c(my.plot.legend.size,
                                     my.plot.legend.face),
                       size = surv.curv.size,
                       censor= my.censor,
                       censor.size=my.censor.size,
                       censor.shape=my.censor.shape,
                       ggtheme = custom_theme())
  
  #Saving the plot to .svg or other specified format
  #with particular parameters
  ggsave(filename = my.km.plot,
         plot     = print(sur_Plot$plot, newpage = FALSE),
         device   = my.km.plot.type, dpi=my.plot.dpi,
         width    = my.plot.width, height = my.plot.height,
         units    = my.plot.unit)
  
  #Returning our finished KM plot----
  return(sur_Plot)
}

