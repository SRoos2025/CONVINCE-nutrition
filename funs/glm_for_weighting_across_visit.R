glm_for_weighting_across_visit <- function(i) {
    dat_imp <- filter(example, visit == i)
    #make model for chance of not being censored, inf_cens == 0
    glm(
        inf_cens == 0 ~ age + sex + dial_vintage + group +
            #time-varying variables
            uf_vol + sbp_pre + urea_pre_mgdl +
            phos_pre_mmoll + residual_urine_out + hdf_convol +
            dbp_pre + ktv_pre + hb_pre_mmoll + crp_pre_mgdl + k_pre_mmoll + SCR_MH_AP_OCCUR:SCR_MH_COPD_OCCUR,
        data   = dat_imp,
        family = binomial() #logistic regression for binary outcome (inf cens 1 or 0)
    )
}

    