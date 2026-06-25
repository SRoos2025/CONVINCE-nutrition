set_endpt_and_npcrscale <- function(data) {
    data <- data %>%
        mutate(
            endpt = if_else(endpt == 2, 0, endpt),
            npcr_scale = npcr *10
        )
}
