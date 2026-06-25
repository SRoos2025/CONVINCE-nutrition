prob_to_ipcw_total <- function(df) {
    map(unique(df[["visit"]]), function(v) {
        df_visit <- filter(df, visit == v)
        
        df_visit %>%
            mutate(
                ps_ipcw = ps,
                ipcw = case_when(
                    inf_cens == 1 ~ 1 / (1 - ps_ipcw),
                    inf_cens == 0 ~ 1 / ps_ipcw
                )
            )
    })
}

