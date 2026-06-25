svy_lme <- function (data, outcome) {
    
    formula <- as.formula(
        paste0(
            outcome,
            " ~ visit_time * group + baseline_", outcome,
            " + (1 | id) + (-1 + visit_time | id)" #for svylme random slope can be independent with -1 + 
        )
    )
    #create dummy variable for which w2 applies
    data$obs <- seq_len(nrow(data)) 
    
    # survey design
    design <- svydesign(
        ids = ~id + obs,
        weights = ~ipcw + w2,
        data = data
    )
    
    # mixed model with survey weights
    fit <- svy2lme(
        formula = formula,
        design = design
    )
    
    coef <- fit$beta
    V <- as.matrix(fit$Vbeta)
    se <- sqrt(diag(vcov(fit)))
    
    return(list(
        model = fit,
        coef = coef,
        SE = se
    ))
}
