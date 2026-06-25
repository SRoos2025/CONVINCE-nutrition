#function to make plots of cox splines of nutr_variables 
nutr_flex_plot <- function(data, nutr_var, x_lab, center_val, break_min = 0, break_max = 3, breaks = 0.5, y_limits = c(0, 5)) {

   plot<-  ggplot(data,
                         aes(x = .data[[nutr_var]],
                             y = hr,
                             ymin = ll,
                             ymax = ul)) +
    geom_hline(yintercept = 1,
               linetype = "dashed",
               colour = "black",
               alpha = 0.4) +
    geom_ribbon(alpha = 0.25,
                fill = "#DC267F") +
    geom_line(size = 1,
              colour = "#DC267F") +
       geom_vline(xintercept = center_val, 
                  linetype = "dashed", 
                  color = "black", 
                  alpha = 0.4) +
    labs(
        x = x_lab,
        y = "Adjusted* Hazard Ratio for mortality"
    ) +
       scale_x_continuous(
           #define second axis
           sec.axis = sec_axis(~.,
                               breaks = center_val,
                               labels = paste0("Reference at ", center_val)),
           #the breaks on x-axis
      breaks = seq(break_min, break_max, breaks)) +
       coord_cartesian(ylim = y_limits) +
    theme_bw(base_size = 20)+
       theme(
                   panel.grid.minor = element_blank(),
                  axis.title.x.top = element_text(colour = "black"),
                   axis.title.y.left = element_text(size = 20),
                  axis.title.x = element_text(size = 20),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_line(color = "grey80"),
                   panel.grid.minor.y = element_line(color = "grey80"),
                   axis.ticks = element_line(color = "black"),
                  axis.text.x.top = element_text(size = 20),
                  axis.text.x.bottom = element_text(size = 20),
                  axis.text.y.left = element_text(size = 20)
               )
return (plot)
}
