#function to pool predictions manually
pooled_pred <- function(data) {
    data <- data %>%
    group_by(data[["x"]], data[["group"]]) %>%
    reframe(
        predicted_mean = mean(predicted),
        within_var = mean(std.error^2, na.rm = TRUE),
        between_var = var(predicted, na.rm = TRUE),
        total_se = sqrt(within_var + (between_var + between_var/m)),
        conf.low = predicted_mean - 1.96 * total_se,
        conf.high = predicted_mean + 1.96 * total_se,
        group = group,
    x =x 
    )%>%
        ungroup()
    return(data)
}
