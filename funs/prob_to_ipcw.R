prob_to_ipcw <- function(i, ps) {
    dat_imp <- filter(data_unweighted, .imp == i)
    
    dat_imp <- dat_imp %>%
        mutate(
            ps_ipcw = ps,
            ipcw = case_when(
                inf_cens == 1 ~ 1 / (1 - ps_ipcw),  
                inf_cens == 0 ~ 1 / ps_ipcw         
            )
        )
    
    dat_imp
}
