pool_manual <- function (data, m) { 
pooled <- data %>%
    group_by(term) %>%
    summarise(
        m = m,
        estimate = mean(estimate, na.rm = TRUE),
        within_var = mean(std.error^2, na.rm = TRUE),
        between_var = var(estimate, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        total_var = within_var + (between_var + between_var/m),
        sd = sqrt(total_var),
        z = estimate / sd,
        conf_low = estimate - 1.96 * sd,
        conf_high = estimate + 1.96 * sd
    ) %>%
    select(term, estimate, sd, conf_low, conf_high, z)

return(pooled)
}
