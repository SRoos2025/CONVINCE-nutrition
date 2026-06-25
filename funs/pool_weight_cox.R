pool_weight_cox <- function(data, m) {
    data %>%
        group_by(term) %>%
        mutate(robust_se = sqrt(diag(vcov(fit))))
        summarise(
            estimate_mean = mean(estimate, na.rm = TRUE),
            
            # binnen-imputatie variantie (robust!)
            within_var = mean(robust_se^2, na.rm = TRUE),
            
            # tussen-imputatie variantie
            between_var = var(estimate, na.rm = TRUE),
            
            # totale variantie volgens Rubin
            total_var = within_var + (1 + 1/m) * between_var,
            total_se = sqrt(total_var),
            
            conf.low = estimate_mean - 1.96 * total_se,
            conf.high = estimate_mean + 1.96 * total_se
        )
}
