#for each outcome, we want to apply the formula visit * group + the baseline value of the outcome, clustered by id
fit_wemix <- function(data, outcome) {
    
    formula <- as.formula(
        paste0(
            outcome,
            " ~ visit_time * group + baseline_", outcome,  #fixed effects
            " +(1 | id)"  #+ (0 + visit_time | id)" #random effects and random slope
        )
    )
    #define weight 1, as WeMix requires 2 weiths minimum
    data[["w2"]] <- 1
    
   mix(
        formula = formula,
        data = data,
        weights = c("ipcw", "w2"))   # only level 1 weights (per individual, no group weights)
}
