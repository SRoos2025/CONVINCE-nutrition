svy_glm <- function (data, outcome) {
    formula <- as.formula(
        paste0(
            outcome,
            " ~ visit_time * group + baseline_", outcome
            #" +(1 | id) + (0 + visit_time | id)" random slope not possible
        )
    )

    
    # define survey design (clustered by id)
    design <- svydesign(
        ids = ~id,
        weights = ~ipcw,
        data = data
    )
    
    fit <- svyglm(
        formula = formula,
        design = design
    )
    #compute standard error
    se <- sqrt(diag(vcov(fit)))
    
    return(list(model = fit, SE = se
    ))
    
}
