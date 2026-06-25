predict_weight <- function(i) {
    dat_imp <- filter(data_unweighted, .imp == i)
    predict(fit_list[[paste0("imp_", i)]],
            #For a binomial model the default predictions are of log-odds (probabilities on logit scale)
            #and type = "response" gives the predicted probabilities.
            type = "response") 
}
