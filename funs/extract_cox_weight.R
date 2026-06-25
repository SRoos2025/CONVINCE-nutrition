extract_cox_weight <- function(fits) {
    bind_rows(
        lapply(fits, function(fit) {
            data.frame(
                term = names(fit$coefficients),
                estimate = fit$coefficients,
                se = sqrt(diag(fit$var))
            )
        }),
        .id = "imp"
    )
}
