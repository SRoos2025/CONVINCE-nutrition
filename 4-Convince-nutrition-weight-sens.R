#nutrition analyses
#weighted, over all imputations
#script by S. Roos 

#0. set-up----
#load packages 
pacman::p_load("conflicted", # Package conflicts 
               "here", # Relative directories 
               "tidyverse", # Data wrangling 
               "patchwork", # Adding plots together
               "magrittr",
               "stats", #for glm model
               "mice", # to pool outcomes using Rubins rules
               "purrr",
               "ggeffects", # to predict values from linear mixed models
               "broom.mixed", # supports mixed models, returns tidy output
               "ggplot2",# to plot the predicted values) # Printing tables to word ) 
               "survival",#to use cox and penalized spline) 
               "gtsummary",#to make a formatted table )
               "WeMix",
               "splines")#for robust standard error with weighted linear mixed model
#solve conflicts 
conflicts_prefer(dplyr::filter) # Between dplyr & stats 
conflicts_prefer(dplyr::select) # Between dplyr & MASS

#set working directory 
path <- "/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/"

setwd("/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/")

#load data. 
#we perform 3 analyses, unweighted, weighted across all imputations and weighted across all imputations AND visits
#across all imputations AND visits (try_out_data_with_ipcw_3 from script 3-Convince_IPCW)
load(paste0(path, "data_weight_allvisits.Rdata"))

# Source relevant functions
walk(list.files(paste0(path, "R/funs")), ~ source(paste0(path, "R/funs/", .x)))

#1.0 datacleaning----
#make visit a factor and #set fixed times for each visit, take the mean of each visit
#you can either use the factor or not in the function, if you want time continuous for the estimates, set to FALSE
data_weight_allvisits <- set_visit_time(data_weight_allvisits, as_factor = TRUE) 
#for all analyses, status can be 2 or 0 there is no seperate competing risk in our analyses 
#define npcr scale, so you can perform linear cox for increase per 0.1 in npcr
#we want the hazard for npcr for each 0.1 increase in npcr, so we multiply exposure by 10 by creating npcrscale
data_weight_allvisits<- set_endpt_and_npcrscale(data_weight_allvisits)

#define baseline also for CRP (and Kt/V for other publication)
data_weight_allvisits <- data_weight_allvisits %>%
    group_by(.imp, id)%>%
    mutate(
        baseline_crp_pre_mgdl = crp_pre_mgdl[visit==0],
        baseline_ktv_pre = ktv_pre[visit == 0]
    )%>%
    ungroup()

#across imputations AND visits
imp_list_weight <- data_weight_allvisits %>%
    group_split(.imp)

##2 LMM nPCR weighted----
#apply LMM function each set in the imputation list
#this is for the weighting fitted on each imputation AND visit
#first we use WeMix
fits_npcr_wemix <- map(imp_list_weight, fit_wemix, outcome = "npcr")

#tidy package does not work with wemix, so we have to create tidy output ourselves using the self written function tidy wemix results
result_npcr_wemix <- tidy_we_mix_results(fits_npcr_wemix)
#pool
pool_npcr_wemix <- rubins_rules_lmm(result_npcr_wemix, m = 10)

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_npcr_wemix <- pool_npcr_wemix %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))

names_matrix <- c(
    "(Intercept)", "Residual", "baseline_npcr", "group1", "id.(Intercept)",
    "visit_time12", "visit_time12:group1",
    "visit_time15", "visit_time15:group1",
    "visit_time18", "visit_time18:group1",
    "visit_time21", "visit_time21:group1",
    "visit_time24", "visit_time24:group1",
    "visit_time27", "visit_time27:group1",
    "visit_time3",
    "visit_time30", "visit_time30:group1",
    "visit_time33", "visit_time33:group1",
    "visit_time36", "visit_time36:group1",
    "visit_time39", "visit_time39:group1",
    "visit_time3:group1",
    "visit_time42", "visit_time42:group1",
    "visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)
#first we set the matrix to 0, we want the number of rows to correspond to plot_data_wemix
matr <- matrix(0, nrow = nrow(plot_data_wemix), ncol = length(names_matrix))
#set names of columns. in matrix to the names of the terms of pool_we_mix
colnames(matr) <- names_matrix

#to be able to check if everything goes well, we give rownames to visit time and group
rownames(matr) <- paste0(
    plot_data_wemix$visit_time,
    plot_data_wemix$group,
    sep = "_"
)

#define mean baseline values of each group
mean_baseline_hd  <- data_weight_allvisits %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_npcr)) %>%
    pull(m)

mean_baseline_hdf <- data_weight_allvisits%>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_npcr)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline npcr value to the corresponding mean value of the group
matr[, "(Intercept)"]    <- 1
matr[, "group1"]         <- plot_data_wemix$group
matr[, "id.(Intercept)"] <- 0

matr[, "baseline_npcr"] <- if_else(plot_data_wemix$group == 0, mean_baseline_hd, mean_baseline_hdf)

#loop over visit times and interaction terms (with :group1)
for (visit_time in unique(plot_data_wemix$visit_time)) {
    
    vt_name <- paste0("visit_time", visit_time)
    vt_int  <- paste0(vt_name, ":group1")
    
    #check for the visit time if it corresponds to the visit time in plot data wemix
    rows <- plot_data_wemix$visit_time == visit_time
    
    #if it corresponds tot the matrix column, set to 1
    if (vt_name %in% colnames(matr)) {
        matr[rows, vt_name] <- 1
    }
    
    #if it matches to the colname of the interaction, and group is 1 then set to 1.
    if (vt_int %in% colnames(matr)) {
        matr[rows & plot_data_wemix$group == 1, vt_int] <- 1
    }
}

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_npcr_wemix$theta
ll <- pool_npcr_wemix$ll
ul <- pool_npcr_wemix$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_npcr_wemix$term
names(ll) <- pool_npcr_wemix$term
names(ul) <- pool_npcr_wemix$term
#just to be sure, we make sure theta ll and ul are in exact same order as the matrix for the multiplication
theta <- theta[colnames(matr)]
ll <- ll[colnames(matr)]
ul <- ul[colnames(matr)]
#matrix multiplication
pred <- matr %*% theta
ll <- matr %*% ll
ul <- matr %*% ul
#assign predictions as a seperate column in plot_data_wemix
plot_data_wemix$pred <- as.numeric(pred)
plot_data_wemix$ll <- as.numeric(ll)
plot_data_wemix$ul <- as.numeric(ul)

plots <- list()
#plot the predicted values            
plots[[1]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "nPCR (g/kg/day)", y_limits = c(0.5, 1.8))
plots[[1]]

##2 LMM BMI weighted----
#apply LMM function each set in the imputation list
#this is for the weighting fitted on each imputation AND visit
#first we use WeMix
fits_bmi_wemix <- map(imp_list_weight, fit_wemix, outcome = "bmi")

#tidy package does not work with wemix, so we have to create tidy output ourselves using the self written function tidy wemix results
result_bmi_wemix <- tidy_we_mix_results(fits_bmi_wemix)
#pool
pool_bmi_wemix <- rubins_rules_lmm(result_bmi_wemix, m = 10)

table_bmi_wemix <- pool_bmi_wemix %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))

names_matrix <- c(
    "(Intercept)", "Residual", "baseline_bmi", "group1", "id.(Intercept)",
    "visit_time12", "visit_time12:group1",
    "visit_time15", "visit_time15:group1",
    "visit_time18", "visit_time18:group1",
    "visit_time21", "visit_time21:group1",
    "visit_time24", "visit_time24:group1",
    "visit_time27", "visit_time27:group1",
    "visit_time3",
    "visit_time30", "visit_time30:group1",
    "visit_time33", "visit_time33:group1",
    "visit_time36", "visit_time36:group1",
    "visit_time39", "visit_time39:group1",
    "visit_time3:group1",
    "visit_time42", "visit_time42:group1",
    "visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)
#first we set the matrix to 0, we want the number of rows to correspond to plot_data_wemix
matr <- matrix(0, nrow = nrow(plot_data_wemix), ncol = length(names_matrix))
#set names of columns. in matrix to the names of the terms of pool_we_mix
colnames(matr) <- names_matrix

#to be able to check if everything goes well, we give rownames to visit time and group
rownames(matr) <- paste0(
    plot_data_wemix$visit_time,
    plot_data_wemix$group,
    sep = "_"
)

#define mean baseline values of each group
mean_baseline_hd  <- data_weight_allvisits %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_bmi)) %>%
    pull(m)

mean_baseline_hdf <- data_weight_allvisits %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_bmi)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline npcr value to the corresponding mean value of the group
matr[, "(Intercept)"]    <- 1
matr[, "group1"]         <- plot_data_wemix$group
matr[, "id.(Intercept)"] <- 0

matr[, "baseline_bmi"] <- if_else(plot_data_wemix$group == 0, mean_baseline_hd, mean_baseline_hdf)

#loop over visit times and interaction terms (with :group1)
for (visit_time in unique(plot_data_wemix$visit_time)) {
    
    vt_name <- paste0("visit_time", visit_time)
    vt_int  <- paste0(vt_name, ":group1")
    
    #check for the visit time if it corresponds to the visit time in plot data wemix
    rows <- plot_data_wemix$visit_time == visit_time
    
    #if it corresponds tot the matrix column, set to 1
    if (vt_name %in% colnames(matr)) {
        matr[rows, vt_name] <- 1
    }
    
    #if it matches to the colname of the interaction, and group is 1 then set to 1.
    if (vt_int %in% colnames(matr)) {
        matr[rows & plot_data_wemix$group == 1, vt_int] <- 1
    }
}

# we can check if everything went well by printing the columnames from the matrix next to the colnames of plot_data_wemix
cbind(
    plot_data_wemix,
    matr[, c("visit_time39", "visit_time39:group1")]
)

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_bmi_wemix$theta
ll <- pool_bmi_wemix$ll
ul <- pool_bmi_wemix$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_bmi_wemix$term
names(ll) <- pool_bmi_wemix$term
names(ul) <- pool_bmi_wemix$term
#just to be sure, we make sure theta ll and ul are in exact same order as the matrix for the multiplication
theta <- theta[colnames(matr)]
ll <- ll[colnames(matr)]
ul <- ul[colnames(matr)]
#matrix multiplication
pred <- matr %*% theta
ll <- matr %*% ll
ul <- matr %*% ul
#assign predictions as a seperate column in plot_data_wemix
plot_data_wemix$pred <- as.numeric(pred)
plot_data_wemix$ll <- as.numeric(ll)
plot_data_wemix$ul <- as.numeric(ul)

expr_bmi <- expression(BMI (kg/m^2))
#plot the predicted values            
plots[[2]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = expr_bmi, y_lim = c(20, 40))
plots[[2]]

##3 LMM LTI weighted----
#apply LMM function each set in the imputation list
#this is for the weighting fitted on each imputation AND visit
#first we use WeMix
fits_lti_wemix <- map(imp_list_weight, fit_wemix, outcome = "lti")

#tidy package does not work with wemix, so we have to create tidy output ourselves using the self written function tidy wemix results
result_lti_wemix <- tidy_we_mix_results(fits_lti_wemix)
#pool
pool_lti_wemix <- rubins_rules_lmm(result_lti_wemix, m = 10)

table_lti_wemix <- pool_lti_wemix %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))

names_matrix <- c(
    "(Intercept)", "Residual", "baseline_lti", "group1", "id.(Intercept)",
    "visit_time12", "visit_time12:group1",
    "visit_time15", "visit_time15:group1",
    "visit_time18", "visit_time18:group1",
    "visit_time21", "visit_time21:group1",
    "visit_time24", "visit_time24:group1",
    "visit_time27", "visit_time27:group1",
    "visit_time3",
    "visit_time30", "visit_time30:group1",
    "visit_time33", "visit_time33:group1",
    "visit_time36", "visit_time36:group1",
    "visit_time39", "visit_time39:group1",
    "visit_time3:group1",
    "visit_time42", "visit_time42:group1",
    "visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)
#first we set the matrix to 0, we want the number of rows to correspond to plot_data_wemix
matr <- matrix(0, nrow = nrow(plot_data_wemix), ncol = length(names_matrix))
#set names of columns. in matrix to the names of the terms of pool_we_mix
colnames(matr) <- names_matrix

#to be able to check if everything goes well, we give rownames to visit time and group
rownames(matr) <- paste0(
    plot_data_wemix$visit_time,
    plot_data_wemix$group,
    sep = "_"
)

#define mean baseline values of each group
mean_baseline_hd  <- data_weight_allvisits %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_lti)) %>%
    pull(m)

mean_baseline_hdf <- data_weight_allvisits %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_lti)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline npcr value to the corresponding mean value of the group
matr[, "(Intercept)"]    <- 1
matr[, "group1"]         <- plot_data_wemix$group
matr[, "id.(Intercept)"] <- 0

matr[, "baseline_lti"] <- if_else(plot_data_wemix$group == 0, mean_baseline_hd, mean_baseline_hdf)

#loop over visit times and interaction terms (with :group1)
for (visit_time in unique(plot_data_wemix$visit_time)) {
    
    vt_name <- paste0("visit_time", visit_time)
    vt_int  <- paste0(vt_name, ":group1")
    
    #check for the visit time if it corresponds to the visit time in plot data wemix
    rows <- plot_data_wemix$visit_time == visit_time
    
    #if it corresponds tot the matrix column, set to 1
    if (vt_name %in% colnames(matr)) {
        matr[rows, vt_name] <- 1
    }
    
    #if it matches to the colname of the interaction, and group is 1 then set to 1.
    if (vt_int %in% colnames(matr)) {
        matr[rows & plot_data_wemix$group == 1, vt_int] <- 1
    }
}

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_lti_wemix$theta
ll <- pool_lti_wemix$ll
ul <- pool_lti_wemix$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_lti_wemix$term
names(ll) <- pool_lti_wemix$term
names(ul) <- pool_lti_wemix$term
#just to be sure, we make sure theta ll and ul are in exact same order as the matrix for the multiplication
theta <- theta[colnames(matr)]
ll <- ll[colnames(matr)]
ul <- ul[colnames(matr)]
#matrix multiplication
pred <- matr %*% theta
ll <- matr %*% ll
ul <- matr %*% ul
#assign predictions as a seperate column in plot_data_wemix
plot_data_wemix$pred <- as.numeric(pred)
plot_data_wemix$ll <- as.numeric(ll)
plot_data_wemix$ul <- as.numeric(ul)

expr_lti <- expression (LTI (mg/m^2))
#plot the predicted values            
plots[[3]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = expr_lti, y_lim = c(15, 22))
plots[[3]]

##4 LMM SCI weighted----
#apply LMM function each set in the imputation list
#this is for the weighting fitted on each imputation AND visit
#first we use WeMix
fits_sci_wemix <- map(imp_list_weight, fit_wemix, outcome = "sci")

#tidy package does not work with wemix, so we have to create tidy output ourselves using the self written function tidy wemix results
result_sci_wemix <- tidy_we_mix_results(fits_sci_wemix)
#pool
pool_sci_wemix <- rubins_rules_lmm(result_sci_wemix , m = 10)

table_sci_wemix <- pool_sci_wemix %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))

names_matrix <- c(
    "(Intercept)", "Residual", "baseline_sci", "group1", "id.(Intercept)",
    "visit_time12", "visit_time12:group1",
    "visit_time15", "visit_time15:group1",
    "visit_time18", "visit_time18:group1",
    "visit_time21", "visit_time21:group1",
    "visit_time24", "visit_time24:group1",
    "visit_time27", "visit_time27:group1",
    "visit_time3",
    "visit_time30", "visit_time30:group1",
    "visit_time33", "visit_time33:group1",
    "visit_time36", "visit_time36:group1",
    "visit_time39", "visit_time39:group1",
    "visit_time3:group1",
    "visit_time42", "visit_time42:group1",
    "visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)
#first we set the matrix to 0, we want the number of rows to correspond to plot_data_wemix
matr <- matrix(0, nrow = nrow(plot_data_wemix), ncol = length(names_matrix))
#set names of columns. in matrix to the names of the terms of pool_we_mix
colnames(matr) <- names_matrix

#to be able to check if everything goes well, we give rownames to visit time and group
rownames(matr) <- paste0(
    plot_data_wemix$visit_time,
    plot_data_wemix$group,
    sep = "_"
)

#define mean baseline values of each group
mean_baseline_hd  <- data_weight_allvisits %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_sci)) %>%
    pull(m)

mean_baseline_hdf <- data_weight_allvisits %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_sci)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline npcr value to the corresponding mean value of the group
matr[, "(Intercept)"]    <- 1
matr[, "group1"]         <- plot_data_wemix$group
matr[, "id.(Intercept)"] <- 0

matr[, "baseline_sci"] <- if_else(plot_data_wemix$group == 0, mean_baseline_hd, mean_baseline_hdf)

#loop over visit times and interaction terms (with :group1)
for (visit_time in unique(plot_data_wemix$visit_time)) {
    
    vt_name <- paste0("visit_time", visit_time)
    vt_int  <- paste0(vt_name, ":group1")
    
    #check for the visit time if it corresponds to the visit time in plot data wemix
    rows <- plot_data_wemix$visit_time == visit_time
    
    #if it corresponds tot the matrix column, set to 1
    if (vt_name %in% colnames(matr)) {
        matr[rows, vt_name] <- 1
    }
    
    #if it matches to the colname of the interaction, and group is 1 then set to 1.
    if (vt_int %in% colnames(matr)) {
        matr[rows & plot_data_wemix$group == 1, vt_int] <- 1
    }
}

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_sci_wemix$theta
ll <- pool_sci_wemix$ll
ul <- pool_sci_wemix$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_sci_wemix$term
names(ll) <- pool_sci_wemix$term
names(ul) <- pool_sci_wemix$term
#just to be sure, we make sure theta ll and ul are in exact same order as the matrix for the multiplication
theta <- theta[colnames(matr)]
ll <- ll[colnames(matr)]
ul <- ul[colnames(matr)]
#matrix multiplication
pred <- matr %*% theta
ll <- matr %*% ll
ul <- matr %*% ul
#assign predictions as a seperate column in plot_data_wemix
plot_data_wemix$pred <- as.numeric(pred)
plot_data_wemix$ll <- as.numeric(ll)
plot_data_wemix$ul <- as.numeric(ul)

#plot the predicted values            
plots[[4]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "SCI (mg/kg/day)", y_lim = c(15, 25))
plots[[4]]

#1.1.1.1 combine LMM plots in one----
figure_1 <- wrap_plots(plots,nrow = 2) +  
    plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect")  # combineer legendes

# marge between plos
figure_1 <- figure_1 & theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2), legend.position = "bottom",  legend.key.size = unit(0.5, "cm"), 
                             legend.text = element_text(size = 20),
                             plot.tag = element_text(size = 30))  
figure_1
#save
ggsave("lmm_weighted_sens_nutr.png", figure_1,
       width = 20, height = 12, dpi = 1200)




##LMM CRP 
fits_crp_wemix <- map(imp_list_weight, fit_wemix, outcome = "crp_pre_mgdl")

result_crp_wemix <- tidy_we_mix_results(fits_crp_wemix)

pool_crp_wemix <- rubins_rules_lmm(result_crp_wemix , m = 10)
selected_terms <- paste0("visit_time", seq(3, 45, by = 3))
table_crp_wemix <- pool_crp_wemix %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))

names_matrix <- c(
    "(Intercept)", "Residual", "baseline_crp_pre_mgdl", "group1", "id.(Intercept)",
    "visit_time12", "visit_time12:group1",
    "visit_time15", "visit_time15:group1",
    "visit_time18", "visit_time18:group1",
    "visit_time21", "visit_time21:group1",
    "visit_time24", "visit_time24:group1",
    "visit_time27", "visit_time27:group1",
    "visit_time3",
    "visit_time30", "visit_time30:group1",
    "visit_time33", "visit_time33:group1",
    "visit_time36", "visit_time36:group1",
    "visit_time39", "visit_time39:group1",
    "visit_time3:group1",
    "visit_time42", "visit_time42:group1",
    "visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)
#first we set the matrix to 0, we want the number of rows to correspond to plot_data_wemix
matr <- matrix(0, nrow = nrow(plot_data_wemix), ncol = length(names_matrix))
#set names of columns. in matrix to the names of the terms of pool_we_mix
colnames(matr) <- names_matrix

#to be able to check if everything goes well, we give rownames to visit time and group
rownames(matr) <- paste0(
    plot_data_wemix$visit_time,
    plot_data_wemix$group,
    sep = "_"
)

#define mean baseline values of each group
mean_baseline_hd  <- data_weight_allvisits %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_crp_pre_mgdl)) %>%
    pull(m)

mean_baseline_hdf <- data_weight_allvisits %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_crp_pre_mgdl)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline npcr value to the corresponding mean value of the group
matr[, "(Intercept)"]    <- 1
matr[, "group1"]         <- plot_data_wemix$group
matr[, "id.(Intercept)"] <- 0

matr[, "baseline_crp_pre_mgdl"] <- if_else(plot_data_wemix$group == 0, mean_baseline_hd, mean_baseline_hdf)

#loop over visit times and interaction terms (with :group1)
for (visit_time in unique(plot_data_wemix$visit_time)) {
    
    vt_name <- paste0("visit_time", visit_time)
    vt_int  <- paste0(vt_name, ":group1")
    
    #check for the visit time if it corresponds to the visit time in plot data wemix
    rows <- plot_data_wemix$visit_time == visit_time
    
    #if it corresponds tot the matrix column, set to 1
    if (vt_name %in% colnames(matr)) {
        matr[rows, vt_name] <- 1
    }
    
    #if it matches to the colname of the interaction, and group is 1 then set to 1.
    if (vt_int %in% colnames(matr)) {
        matr[rows & plot_data_wemix$group == 1, vt_int] <- 1
    }
}

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_crp_wemix$theta
ll <- pool_crp_wemix$ll
ul <- pool_crp_wemix$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_crp_wemix$term
names(ll) <- pool_crp_wemix$term
names(ul) <- pool_crp_wemix$term
#just to be sure, we make sure theta ll and ul are in exact same order as the matrix for the multiplication
theta <- theta[colnames(matr)]
ll <- ll[colnames(matr)]
ul <- ul[colnames(matr)]
#matrix multiplication
pred <- matr %*% theta
ll <- matr %*% ll
ul <- matr %*% ul
#assign predictions as a seperate column in plot_data_wemix
plot_data_wemix$pred <- as.numeric(pred)
plot_data_wemix$ll <- as.numeric(ll)
plot_data_wemix$ul <- as.numeric(ul)

#lower limit can not be below 0
plot_data_wemix <- plot_data_wemix %>%
    mutate(
        ll = if_else(ll<0, 0, ll)
    )

crp_plot <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "CRP (mg/dL)", y_lim = c(0, 10))
crp_plot

ggsave("crp_plot_weight_sens.png", crp_plot,
       width = 20, height = 12, dpi = 600)

#LM Kt/V----
fits_ktv_wemix <- map(imp_list_weight, fit_wemix, outcome = "ktv_pre")

result_ktv_wemix <- tidy_we_mix_results(fits_ktv_wemix)

pool_ktv_wemix <- rubins_rules_lmm(result_ktv_wemix , m = 10)

table_ktv_wemix <- pool_ktv_wemix %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))



names_matrix <- c(
    "(Intercept)", "Residual", "baseline_ktv_pre", "group1", "id.(Intercept)",
    "visit_time12", "visit_time12:group1",
    "visit_time15", "visit_time15:group1",
    "visit_time18", "visit_time18:group1",
    "visit_time21", "visit_time21:group1",
    "visit_time24", "visit_time24:group1",
    "visit_time27", "visit_time27:group1",
    "visit_time3",
    "visit_time30", "visit_time30:group1",
    "visit_time33", "visit_time33:group1",
    "visit_time36", "visit_time36:group1",
    "visit_time39", "visit_time39:group1",
    "visit_time3:group1",
    "visit_time42", "visit_time42:group1",
    "visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)
#first we set the matrix to 0, we want the number of rows to correspond to plot_data_wemix
matr <- matrix(0, nrow = nrow(plot_data_wemix), ncol = length(names_matrix))
#set names of columns. in matrix to the names of the terms of pool_we_mix
colnames(matr) <- names_matrix

#to be able to check if everything goes well, we give rownames to visit time and group
rownames(matr) <- paste0(
    plot_data_wemix$visit_time,
    plot_data_wemix$group,
    sep = "_"
)

#define mean baseline values of each group
mean_baseline_hd  <- data_weight_allvisits %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_ktv_pre)) %>%
    pull(m)

mean_baseline_hdf <- data_weight_allvisits %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_ktv_pre)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline npcr value to the corresponding mean value of the group
matr[, "(Intercept)"]    <- 1
matr[, "group1"]         <- plot_data_wemix$group
matr[, "id.(Intercept)"] <- 0

matr[, "baseline_ktv_pre"] <- if_else(plot_data_wemix$group == 0, mean_baseline_hd, mean_baseline_hdf)

#loop over visit times and interaction terms (with :group1)
for (visit_time in unique(plot_data_wemix$visit_time)) {
    
    vt_name <- paste0("visit_time", visit_time)
    vt_int  <- paste0(vt_name, ":group1")
    
    #check for the visit time if it corresponds to the visit time in plot data wemix
    rows <- plot_data_wemix$visit_time == visit_time
    
    #if it corresponds tot the matrix column, set to 1
    if (vt_name %in% colnames(matr)) {
        matr[rows, vt_name] <- 1
    }
    
    #if it matches to the colname of the interaction, and group is 1 then set to 1.
    if (vt_int %in% colnames(matr)) {
        matr[rows & plot_data_wemix$group == 1, vt_int] <- 1
    }
}

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_ktv_wemix$theta
ll <- pool_ktv_wemix$ll
ul <- pool_ktv_wemix$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_ktv_wemix$term
names(ll) <- pool_ktv_wemix$term
names(ul) <- pool_ktv_wemix$term
#just to be sure, we make sure theta ll and ul are in exact same order as the matrix for the multiplication
theta <- theta[colnames(matr)]
ll <- ll[colnames(matr)]
ul <- ul[colnames(matr)]
#matrix multiplication
pred <- matr %*% theta
ll <- matr %*% ll
ul <- matr %*% ul
#assign predictions as a seperate column in plot_data_wemix
plot_data_wemix$pred <- as.numeric(pred)
plot_data_wemix$ll <- as.numeric(ll)
plot_data_wemix$ul <- as.numeric(ul)

ktv_plot <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "Kt/V", y_lim = c(1, 2.5))
ktv_plot

ggsave("ktv_weight_sens.png", ktv_plot,
       width = 16, height = 10, dpi = 600)


#5.0 linear cox weighted over all imputations----
#5.1 npcr----
fits_npcr_cox_weight <- map(imp_list_weight, \(x)cox_nutr(x, "npcr_scale"))
#extract results and put in dataframe
result_cox_weight <- extract_cox_weight(fits_npcr_cox_weight)
#manually pool using rubins rules
#the function for pooling for lmm works the same for the cox
pooled_result_cox_weight <- rubins_rules_linear_cox(result_cox_weight, m = 10)

#5.2 bmi----
#linear cox weighted over all imputations
fits_bmi_cox_weight <- map(imp_list_weight, \(x)cox_nutr(x, "bmi"))
#extract results and put in dataframe
result_cox_bmi_weight <- extract_cox_weight(fits_bmi_cox_weight)
#manually pool using rubins rules
#the function for pooling for lmm works the same for the cox
pooled_result_bmi_cox_weight <- rubins_rules_linear_cox(result_cox_bmi_weight, m = 10)

#5.3 sci----
#linear cox weighted over all imputations
fits_sci_cox_weight <- map(imp_list_weight, \(x)cox_nutr(x, "sci"))
#extract results and put in dataframe
result_cox_sci_weight <- extract_cox_weight(fits_sci_cox_weight)
#manually pool using rubins rules
#the function for pooling for lmm works the same for the cox
pooled_result_sci_cox_weight <- rubins_rules_linear_cox(result_cox_sci_weight, m = 10)

#5.4 lti----
#linear cox weighted over all imputations
fits_lti_cox_weight <- map(imp_list_weight, \(x)cox_nutr(x, "lti"))
#extract results and put in dataframe
result_cox_lti_weight <- extract_cox_weight(fits_lti_cox_weight)
#manually pool using rubins rules
#the function for pooling for lmm works the same for the cox
pooled_result_lti_cox_weight <- rubins_rules_linear_cox(result_cox_lti_weight, m = 10)

#6.0 cox splines----
#6.1 npcr----
fit_spline_npcr_cox_weight <- map(imp_list_weight, \(x)cox_spline_nutr(x, nutr_var = "npcr", 
                                                                           lower_bound = min(data_weight_allvisits[["npcr"]]),
                                                                           upper_bound = max(data_weight_allvisits[["npcr"]]), weights = TRUE)
)
#turn it from list to dataframe
fit_spline_npcr_cox <- bind_rows(fit_spline_npcr_cox_weight)
#apply rubin's rules
fit_spline_npcr_cox <- rubin_rule_cox_spline(fit_spline_npcr_cox, nutr_var = "npcr", center_at = 1.1)
#make a list for the plots
plots_spline_cox <- list()
#plot it
plots_spline_cox[[1]]<- nutr_flex_plot(fit_spline_npcr_cox, nutr_var = "npcr", x_lab = "nPCR (g/kg/day)", center_val = 1.1)

plots_spline_cox[[1]]

#6.2 BMI ----
#weighted per visit AND imputation
fit_spline_bmi_cox_weight <- map(imp_list_weight, \(x)cox_spline_nutr(x, nutr_var = "bmi", 
                                                                          lower_bound = min(data_weight_allvisits[["bmi"]]),
                                                                          upper_bound = 50, weights = TRUE)
)
#turn it from list to dataframe
fit_spline_bmi_cox <- bind_rows(fit_spline_bmi_cox_weight)
#apply rubin's rules
fit_spline_bmi_cox <- rubin_rule_cox_spline(fit_spline_bmi_cox, nutr_var = "bmi", center_at = 25)
expr_bmi <- expression(BMI (kg/m^2))
#plot it
plots_spline_cox[[2]]<- nutr_flex_plot(fit_spline_bmi_cox, nutr_var = "bmi", x_lab = expr_bmi, center_val = 25, break_min = 10, break_max = 50, breaks = 10)

plots_spline_cox[[2]]

#6.3 LTI ----
#weighted per visit AND imputation
#reference value is median of total group
fit_spline_lti_cox_weight <- map(imp_list_weight, \(x)cox_spline_nutr(x, nutr_var = "lti", 
                                                                          lower_bound = min(data_weight_allvisits[["lti"]]),
                                                                          upper_bound = max(data_weight_allvisits[["lti"]]), weights = TRUE)
)
#turn it from list to dataframe
fit_spline_lti_cox <- bind_rows(fit_spline_lti_cox_weight)
#apply rubin's rules
fit_spline_lti_cox <- rubin_rule_cox_spline(fit_spline_lti_cox, nutr_var = "lti", center_at = 17.2)
expr_lti <- expression (LTI (mg/m^2))
#plot it
plots_spline_cox[[3]]<- nutr_flex_plot(fit_spline_lti_cox, nutr_var = "lti", x_lab = expr_lti, center_val = 17.2, break_min = 10, break_max = 40, breaks = 10)

plots_spline_cox[[3]]


#6.3 SCI ----
#weighted per visit AND imputation
fit_spline_sci_cox_weight <- map(imp_list_weight, \(x)cox_spline_nutr(x, nutr_var = "sci", 
                                                                          lower_bound = min(data_weight_allvisits[["sci"]]),
                                                                          upper_bound = max(data_weight_allvisits[["sci"]]), weights = TRUE)
)
#turn it from list to dataframe
fit_spline_sci_cox <- bind_rows(fit_spline_sci_cox_weight)
#apply rubin's rules
fit_spline_sci_cox <- rubin_rule_cox_spline(fit_spline_sci_cox, nutr_var = "sci", center_at = 19.3)
#plot it
plots_spline_cox[[4]]<- nutr_flex_plot(fit_spline_sci_cox, nutr_var = "sci", x_lab = "SCI (mg/kg/day)", center_val = 19.3, break_min = 15, break_max = 30, breaks = 5)

plots_spline_cox[[4]]


figure_1 <- wrap_plots(plots_spline_cox) + plot_annotation(tag_levels = "A")
#check
figure_1 
#save
ggsave("flex_spline_nutr_weighted_sens.png", figure_1, width= 14, height = 10, dpi = 1200)
