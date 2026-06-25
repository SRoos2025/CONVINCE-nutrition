#make distribution plot of weighting for specific study visit
plot_ipcw_visit_unweighted <- function(data) {
   
plot <- ggplot(data, aes(x = ps_ipcw,  fill = as.factor(inf_cens))) +
    #alpha is transparency, 0 is transparent 1 = opaque
    #color NA we dont want lines around the curve
    geom_density(alpha = 0.3, color = NA) +
    facet_wrap(~ visit, ncol = 2, scales = "free_y")+ #max 2 plots per row labeller = labeller(visit = function(x) paste("Visit", x))) + 
    #coord_cartesian(xlim = c(0.8, 1.0)) + # zoom
    labs( x = "Propensity score", y = "Density", title = "Unweighted Distribution of propensity scores per study visit", fill = "Censoring status"
         )+
    theme_minimal()+
    theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    scale_fill_manual(
        values = c("0" = "#648FFF", "1" ="#FFB000"),
        labels = c("0" = "No informative censoring", "1"= "Informative censoring"))
return(plot)
}

