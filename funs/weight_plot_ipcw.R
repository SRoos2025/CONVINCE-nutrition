#make a density plot comparing the propensity score of censored with uncensored with weighting, so weight = ipcw
weight_plot_ipcw <- function(data) {
    ggplot(data, aes(x = ps_ipcw, weight = ipcw, fill = as.factor(inf_cens))) +
        #alpha is transparency, 0 is transparent 1 = opaque
        #color NA we dont want lines around the curve
        geom_density(alpha = 0.3, color = NA) +
        facet_wrap(~ .imp, ncol = 2, #max 2 plots per row
                   labeller = labeller(.imp = function(x) paste("Imputation Iteration", x))) + 
        coord_cartesian(xlim = c(0.6, 1.0)) + # zoom
        labs( x = "Propensity score", y = "Density", title = "Weighted Distribution of propensity scores per imputation", fill = "Censoring status") +
        theme_minimal()+
        theme(
            axis.text.y = element_text(size = 4)) +  
        scale_fill_manual(
            values = c("0" = "#648FFF", "1" ="#FFB000"),
            labels = c("0" = "No informative censoring", "1"= "Informative censoring"))
}
