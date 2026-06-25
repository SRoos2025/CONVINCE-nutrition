cox_spline_mgus <- function (data, lower_bound, upper_bound) {
    #left hand side of formula
    lhs <- "Surv(futime, death)~"
    rhs <- paste0("pspline(age, df = 4)",
                  #degrees of freedom controls how flexible the curve can be. We choose 4 which is most common
                  "+ sex")
    #fit the cox formula
    fit <- coxph(as.formula(paste0(lhs, rhs)), data = data)
    
    #define lower and upper bound
    # lower_bound <- min(data[[nutr_var]])
    # upper_bound <- max(data[[nutr_var]])
    
    # Get grid of observed values, with centred values for non x-axis variables
    dat_grid <- data.frame(
        x = seq(lower_bound, upper_bound, length.out = 500)
    )
    
    names(dat_grid)[1] <- "age"
    
    dat_grid$sex <- "female"
    #get predicitons for fit
    lst_predictions <- predict(fit, newdata = dat_grid, se.fit = TRUE)
    
    # Add predictions to grid
    dat_grid %<>% mutate(
        b  = lst_predictions[[1]],
        se = lst_predictions[[2]]
    )
    
    return(dat_grid)
}
