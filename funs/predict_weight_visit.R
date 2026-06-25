predict_weight_visit <- function(i) {
    dat_imp <- filter(example, visit == i)
    predict(fit_list2[[paste0("visit_", i)]],
            #For a binomial model the default predictions are of log-odds (probabilities on logit scale)
            #and type = "response" gives the predicted probabilities.
            type = "response") 
}
