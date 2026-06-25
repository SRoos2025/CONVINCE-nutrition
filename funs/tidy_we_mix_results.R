#we mix is not supported with tidy
#this function extracts the result of interest from the model fits
#so later on they can be pooled using rubin's rules
#this function is copied and altered from https://raw.githubusercontent.com/flh3/pubdata/main/mixPV/wemix_modelsummary.R
tidy_we_mix_results <- function(fit_list) {
    all_results <- lapply(seq_along(fit_list), function(i) {
        m1 <- fit_list[[i]]
        
        
        re.tmp <- m1$varDF # extract random effect ifnormation
        estimates <- unname(c(m1$coef, re.tmp$vcov)) #fixed + random estimates
        ses <- unname(c(m1$SE, re.tmp$SEvcov)) #corresponding SE's
        
        #create dataframe per fit
        df <- data.frame(
            imp = i, 
            term = c(names(m1$coef), re.tmp$fullGroup),
            estimate = estimates,
            se = ses
        )
        return(df)
    })
    
    #combine into one dataframe
    final_df <- do.call(rbind, all_results)
    return(final_df)
}
