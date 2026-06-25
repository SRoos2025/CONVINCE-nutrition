#make cox penalized spline function for nutritional variables
cox_spline_nutr <- function (data, nutr_var, lower_bound, upper_bound, weights = FALSE) {
    #left hand side of formula
    lhs <- "Surv(t_start, t_stop, endpt)~"
    rhs <- paste0("ns(", nutr_var, ", df = 3)",
#degrees of freedom controls how flexible the curve can be. We choose 4 which is most common
                  "+ age +sex +dial_vintage +phos_pre_mmoll +hb_pre_mmoll +crp_pre_mgdl +pth_pre_pmoll+na_pre_mmoll+k_pre_mmoll +mg_pre_mmoll +SCR_MH_COPD_OCCUR +SCR_MH_DIAB_OCCUR +SCR_MH_CHF_OCCUR")
   #fit the cox formula
    if(weights) {
     fit <- coxph(as.formula(paste0(lhs, rhs)), data = data, 
                  robust = TRUE, 
                  weights = ipcw, 
                  cluster = id)
     } else {
         fit <- coxph(as.formula(paste0(lhs, rhs)), data = data, 
                      cluster = id)
     }
# Get grid of observed values, with centred values for non x-axis variables
    dat_grid <- data.frame(
        x = seq(lower_bound, upper_bound, length.out = 500)
    )
    
    names(dat_grid)[1] <- nutr_var
    
    dat_grid$age <- mean(data[["age"]], na.rm = TRUE)
    dat_grid$sex <- 0
    dat_grid$dial_vintage <- mean(data[["dial_vintage"]], na.rm = TRUE)
    dat_grid$phos_pre_mmoll <- mean(data[["phos_pre_mmoll"]], na.rm = TRUE)
    dat_grid$hb_pre_mmoll <- mean(data[["hb_pre_mmoll"]], na.rm = TRUE)
    dat_grid$crp_pre_mgdl <- mean(data[["crp_pre_mgdl"]], na.rm = TRUE)
    dat_grid$pth_pre_pmoll <- mean(data[["pth_pre_pmoll"]], na.rm = TRUE)
    dat_grid$na_pre_mmoll <- mean(data[["na_pre_mmoll"]], na.rm = TRUE)
    dat_grid$k_pre_mmoll <- mean(data[["k_pre_mmoll"]], na.rm = TRUE)
    dat_grid$mg_pre_mmoll <- mean(data[["mg_pre_mmoll"]], na.rm = TRUE)
    dat_grid$SCR_MH_COPD_OCCUR <- 0
    dat_grid$SCR_MH_DIAB_OCCUR <- 0
    dat_grid$SCR_MH_CHF_OCCUR <- 0
    #get predicitons for fit
    lst_predictions <- predict(fit, newdata = dat_grid, se.fit = TRUE)
    
    # Add predictions to grid
    dat_grid %<>% mutate(
        b  = lst_predictions[[1]],
        se = lst_predictions[[2]]
    )
    
    return(dat_grid)
}
