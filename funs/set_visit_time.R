set_visit_time<- function(data, as_factor) {
    data <- data %>%
        mutate(
        visit_time = case_when(
            visit == 0 ~ 0,
            visit == 1 ~ 3,
            visit == 2 ~ 6,
            visit == 3 ~ 9,
            visit == 4 ~ 12,
            visit == 5 ~ 15,
            visit == 6 ~ 18,
            visit == 7 ~ 21,
            visit == 8 ~ 24,
            visit == 9 ~ 27,
            visit == 10 ~ 30,
            visit == 11 ~ 33,
            visit == 12 ~ 36,
            visit == 13 ~ 39,
            visit == 14 ~ 42,
            visit == 15 ~ 45),
        visit_time = if(as_factor) as.factor(visit_time) else visit_time #make a factor to check non-linear associations
    ) %>%
    ungroup()
}
