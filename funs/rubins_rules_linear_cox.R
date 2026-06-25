#for dataframe with estimate, imp, standard error and estimates
#define length of imputations when calling the formula
rubins_rules_linear_cox <- function(data, m) {
    
    pooled <- data %>%
        group_by(term) %>%
        summarise(
            theta = mean(estimate),                      # mean beta/estimate
            vw = mean(se^2),                  # within imputation variance
            vb =  var(estimate)  # between-imputation variance
        ) %>%
        ungroup() %>%
        mutate(
            se_total = sqrt(vw + vb + vb / m),
            hr = exp(theta),
            ll = exp(theta - 1.96 * se_total),
            ul = exp(theta + 1.96 * se_total)
        )
}
