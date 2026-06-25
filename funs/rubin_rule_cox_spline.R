
#function for Rubin's rules
rubin_rule_cox_spline <- function(data, nutr_var, center_at) {
    
    #with approx function, you return a list of points which linearly interpolate given data points
    #x and y are the numeric vectorns giving the coordinaties of the points to be interpolated
    #xout is an optional set of numeric values specifying where interpolation is to take place
    
            center <- approx(x = data[[nutr_var]], y = data$b, xout = center_at)$y



    data <- data %>%
    arrange(.data[[nutr_var]]) %>% #npcr
    group_by(.data[[nutr_var]]) %>% #npcr
    mutate(
        theta = mean(b - center),
        vw = mean(se^2),
        vb = var(b-center)
    ) %>%
    slice(1L) %>%
    ungroup() %>%
    mutate(
        se_total = sqrt(vw + vb + vb/10),
        hr = exp(theta),
        ll = exp(theta - 1.96 * se_total),
        ul = exp(theta + 1.96 * se_total)
    )
    
}
