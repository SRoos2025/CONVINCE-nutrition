predict_weight_total <- function(i) {
    dat_imp <- filter(data_unweighted, .imp == i)
    
    # map over visits and put in dataframe using bind rows
     bind_rows(
        map(visits, function(v) {
            dat_imp_visit <- filter(dat_imp, visit == v)
            
            dat_imp_visit <- dat_imp_visit %>%
                mutate(
                .imp = i,
                visit = v,
                ps = predict(
                    fit_list3[[paste0("imp_", i)]][[as.character(v)]],
                    type = "response"
                )
            )
        })
    )
}
