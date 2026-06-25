#run cox on nutritional variables
cox_nutr <- function(.data, nutr_var, robust = TRUE, weights = TRUE) {
    lhs <- "Surv(t_start, t_stop, endpt) ~" 
    rhs <-    paste0(nutr_var, 
              #Adjustment is necessary as nutrition is not a randomised intervention
              "+ age +sex +dial_vintage +phos_pre_mmoll +hb_pre_mmoll +crp_pre_mgdl +pth_pre_pmoll+na_pre_mmoll+k_pre_mmoll +mg_pre_mmoll +SCR_MH_COPD_OCCUR +SCR_MH_DIAB_OCCUR +SCR_MH_CHF_OCCUR")
fit <- coxph(as.formula(paste0(lhs, rhs)), data = .data, 
             robust = robust, 
             weights = if(weights) { 
                 ipcw
             } else if (weights == FALSE) {
                 NULL
             }else {
                 NULL
             },
             cluster = id)

return(fit)
    }




