plot_lmm_nutrition <- function(data, y_lab, y_limits) {
    ggplot(data, aes(x = x, y = predicted_mean, color = factor(group), group = factor(group))) +
    # line per group including size
    geom_line(size = 1) +
    #measurementpoints
    geom_point(size = 2) +
    #shadow for confidence interval, alpha 0.2 makes it semi transparent
    #color = NA no border around confidence interval
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(group)), #achter andere aes zetten
                alpha = 0.2, color = NA) +
    labs(
        x = "Time (months)",
        y = y_lab,
        #no text like "group" above the legend
        color = NULL,
        fill = NULL
    ) +
    scale_color_manual(
        values = c("0" = "#648FFF", "1" = "#FE6100"),
        labels = c("0" = "High-flux hemodialysis", "1" = "High-dose hemodiafiltration")
    ) +
    scale_fill_manual(
        values = c("0" = "#648FFF", "1" = "#FE6100"),
        labels = c("0" = "High-flux hemodialysis", "1" = "High-dose hemodiafiltration")
    ) +
    #use coord cartesian so that you only plot within limits but do not throw any datapoints away
    coord_cartesian(ylim = y_limits)+ #
    theme_minimal(base_size = 20) +
    theme(
        panel.border = element_rect(color = "black", size = 1), #black line surrounding plot
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        panel.grid.major.x = element_blank(),  # no vertical lines
        panel.grid.minor.x = element_blank(),  
        panel.grid.major.y = element_line(color = "grey80"),  # horiztontal lines
        panel.grid.minor.y = element_line(color = "grey80"), #fine horizontal lines inbetween 
        axis.ticks = element_line(color = "black"),
        #plot.background = element_blank(),
        axis.title.y.left = element_text(size = 30),
        axis.title.x.bottom = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20 )
       
    )
}
