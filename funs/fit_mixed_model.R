#for each outcome, we want to apply the formula visit * group + the baseline value of the outcome, clustered by id
fit_mixed_model <- function(data, outcome) {
    formula <- as.formula(
        paste0(outcome, " ~ visit_time * group + baseline_", outcome, " + (1|id)")
    )
    lmer(formula, data = data)
}
