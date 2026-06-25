#subgroup analyses for poor nutritional status
#script by S. Roos 
#last updated 12-06-2026

#A weighted over all impuations AND visits
#0. set-up----
#load packages 
pacman::p_load("conflicted", # Package conflicts 
               "here", # Relative directories 
               "tidyverse", # Data wrangling 
               "patchwork", # Adding plots together
               "magrittr", #data cleaning
               "stats", #for glm model
               "mice", # to pool outcomes using Rubins rules
               "purrr", #data cleaning
               "ggeffects", # to predict values from linear mixed models
               "broom.mixed", # supports mixed models, returns tidy output
               "ggplot2",# to plot the predicted values) # Printing tables to word ) 
               "survival",#to use cox and penalized spline) 
               "gtsummary",#to make a formatted table )
               "WeMix", #for weighted linear mixed models
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
load(paste0(path, "data_weight_allvisits_imp.Rdata"))
# Source relevant functions
walk(list.files(paste0(path, "R/funs")), ~ source(paste0(path, "R/funs/", .x)))

#1.0 datacleaning----

#we did not calculate baseline_crp yet
data_weight_allvisits_imp <- data_weight_allvisits_imp %>%
    group_by(.imp, id)%>%
    mutate(
        baseline_crp_pre_mgdl = crp_pre_mgdl[visit==0],
        baseline_ktv_pre = ktv_pre[visit == 0]
    )%>%
    ungroup()

#make visit a factor and #set fixed times for each visit, take the mean of each visit
#you can either use the factor or not in the function, if you want time continuous for the estimates, set to FALSE
data_weight_allvisits_imp <- set_visit_time(data_weight_allvisits_imp, as_factor = TRUE) 
#for all analyses, status can be 2 or 0 there is no seperate competing risk in our analyses 
#define npcr scale, so you can perform linear cox for increase per 0.1 in npcr
#we want the hazard for npcr for each 0.1 increase in npcr, so we multiply exposure by 10 by creating npcrscale
data_weight_allvisits_imp<- set_endpt_and_npcrscale(data_weight_allvisits_imp)
#for subgroup analyses, derive id's for patients with low nutritional status in unimputed data
#load imputed data
load(paste0(path, "imputed_data_convince.Rdata"))
#filter imp 0
unimputed <- dat_imputed %>%
    filter(.imp ==0)

unimputed <- unimputed %>%
    group_by(id)%>%
    mutate(
        baseline_sci = if_else(visit == 0, sci, NA),
        baseline_lti = if_else(visit == 0, lti, NA)
    )%>%
    ungroup()
#for SCI and LTI, we take values below 40th percentile
quantile(unimputed[["baseline_sci"]], probs = 0.40, na.rm = TRUE) # 18.9
quantile(unimputed[["baseline_lti"]], probs = 0.40, na.rm = TRUE) #17.1

unimputed <- unimputed %>%
    group_by(id) %>%
    mutate(
        #define poor values at baseline
        poor_npcr = case_when(
            npcr <0.8 & visit == 0 ~ 1,
            .default = 0),
        poor_npcr = max(poor_npcr, na.rm = TRUE),
        #because all NA is possible, this may return infinite, which we return to 0
        poor_npcr =if_else(is.infinite(poor_npcr), 0, poor_npcr),
        poor_bmi = case_when(
            bmi < 20 & visit == 0 ~ 1,
            .default = 0),
        poor_bmi = max(poor_bmi, na.rm = TRUE),
        poor_bmi =if_else(is.infinite(poor_bmi), 0, poor_bmi),
        poor_lti = case_when(
            lti < 17.1 & visit == 0 ~ 1,
            .default = 0),
        poor_lti = max(poor_lti, na.rm = TRUE),
        poor_lti =if_else(is.infinite(poor_lti), 0, poor_lti),
        poor_sci = case_when(
            sci < 18.9 & visit == 0 ~ 1,
            .default = 0),
        poor_sci = max(poor_sci, na.rm = TRUE),
        poor_sci =if_else(is.infinite(poor_sci), 0, poor_sci))%>%
    ungroup()

poor_nutr_npcr <- unimputed %>%
    filter(poor_npcr==1)
n_distinct(unique(poor_nutr_npcr$id))#148

npcr_ids <- poor_nutr_npcr %>%
    distinct(id) %>%
    pull(id)

poor_nutr_bmi <- unimputed %>%
    filter(poor_bmi==1)
n_distinct(unique(poor_nutr_bmi$id)) #76

bmi_ids <- poor_nutr_bmi %>%
    distinct(id) %>%
    pull(id)

poor_nutr_lti <- unimputed %>%
    filter(poor_lti==1)
n_distinct(unique(poor_nutr_lti$id)) #452

lti_ids <- poor_nutr_lti %>%
    distinct(id) %>%
    pull(id)

poor_nutr_sci <- unimputed %>%
    filter(poor_sci == 1)
n_distinct(unique(poor_nutr_sci$id)) #453

sci_ids <- poor_nutr_sci %>%
    distinct(id) %>%
    pull(id)

#filter these ID's from dat_weight_all visits
poor_npcr <- data_weight_allvisits_imp %>%
    filter(id %in% npcr_ids)

poor_bmi <- data_weight_allvisits_imp %>%
    filter(id %in% bmi_ids)

poor_lti <- data_weight_allvisits_imp %>%
    filter(id %in% lti_ids)

poor_sci <- data_weight_allvisits_imp %>%
    filter(id %in% sci_ids)

#save inbetween         
save(poor_npcr, file = paste0(path, "poor_npcr.Rdata"))
save(poor_bmi, file = paste0(path, "poor_bmi.Rdata"))
save(poor_lti, file = paste0(path, "poor_lti.Rdata"))
save(poor_sci, file = paste0(path, "poor_sci.Rdata"))

imp_list_npcr <- poor_npcr %>%
    group_split(.imp)

imp_list_bmi <- poor_bmi %>%
    group_split(.imp)

imp_list_sci <- poor_sci %>%
    group_split(.imp)

imp_list_lti<- poor_lti %>%
    group_split(.imp)

#A weighted over all imputations AND visits
#2.0 SUBGROUP LMM npcr----
fits_poor_npcr <- map(imp_list_npcr, fit_wemix, outcome = "npcr")

result_poor_npcr <- tidy_we_mix_results(fits_poor_npcr)
#pool
pool_poor_npcr <- rubins_rules_lmm(result_poor_npcr, m = 10)

#remove visit 45 as there is only 1 group left at that point
pool_poor_npcr <- pool_poor_npcr %>%
    filter(term != "visit_time45")

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_npcr <- pool_poor_npcr %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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
    #remove visit 45 as there is only 1 group present at this point
    #"visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_npcr %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_npcr)) %>%
    pull(m)

mean_baseline_hdf <- poor_npcr %>%
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
theta <- pool_poor_npcr$theta
ll <- pool_poor_npcr$ll
ul <- pool_poor_npcr$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_npcr$term
names(ll) <- pool_poor_npcr$term
names(ul) <- pool_poor_npcr$term
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

plots_sub <- list()
#plot the predicted values            
plots_sub[[1]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "nPCR (g/kg/day)", y_limits = c(0.3, 1.6))
plots_sub[[1]]


#3.0 SUBGROUP LMM BMI----
fits_poor_bmi <- map(imp_list_bmi, fit_wemix, outcome = "bmi")

result_poor_bmi <- tidy_we_mix_results(fits_poor_bmi)
#pool
pool_poor_bmi <- rubins_rules_lmm(result_poor_bmi, m = 10)
#note there is no viist 45 available, no persons in this subgroup left over at that timepoint

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_bmi <- pool_poor_bmi %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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
    #remove visit 45
    #"visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_bmi %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_bmi)) %>%
    pull(m)

mean_baseline_hdf <- poor_bmi %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_bmi)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline bmi value to the corresponding mean value of the group
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

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_poor_bmi$theta
ll <- pool_poor_bmi$ll
ul <- pool_poor_bmi$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_bmi$term
names(ll) <- pool_poor_bmi$term
names(ul) <- pool_poor_bmi$term
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
#lower then 10 not possible
plot_data_wemix$ll <- if_else(plot_data_wemix$ll <10, 10, plot_data_wemix$ll)

expr_bmi <- expression(BMI (kg/m^2))
#plot the predicted values            
plots_sub[[2]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = expr_bmi, y_limits = c(10, 40))
plots_sub[[2]]

#4.0 SUBGROUP LMM LTI----
fits_poor_lti <- map(imp_list_lti, fit_wemix, outcome = "lti")

result_poor_lti <- tidy_we_mix_results(fits_poor_lti)
#pool
pool_poor_lti <- rubins_rules_lmm(result_poor_lti, m = 10)

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_lti <- pool_poor_lti %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_lti %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_lti)) %>%
    pull(m)

mean_baseline_hdf <- poor_lti %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_lti)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline lti value to the corresponding mean value of the group
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
theta <- pool_poor_lti$theta
ll <- pool_poor_lti$ll
ul <- pool_poor_lti$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_lti$term
names(ll) <- pool_poor_lti$term
names(ul) <- pool_poor_lti$term
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
plots_sub[[3]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = expr_lti, y_limits = c(10, 22))
plots_sub[[3]]


#5.0 SUBGROUP LMM SCI----
fits_poor_sci <- map(imp_list_sci, fit_wemix, outcome = "sci")

result_poor_sci <- tidy_we_mix_results(fits_poor_sci)
#pool
pool_poor_sci <- rubins_rules_lmm(result_poor_sci, m = 10)

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_sci <- pool_poor_sci %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_sci %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_sci)) %>%
    pull(m)

mean_baseline_hdf <- poor_sci %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_sci)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline sci value to the corresponding mean value of the group
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
theta <- pool_poor_sci$theta
ll <- pool_poor_sci$ll
ul <- pool_poor_sci$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_sci$term
names(ll) <- pool_poor_sci$term
names(ul) <- pool_poor_sci$term
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
plots_sub[[4]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "SCI (mg/kg/day)", y_limits = c(10, 25))
plots_sub[[4]]

#COMBINE SUBGROUPS in 1 plot----
figure_2 <- wrap_plots(plots_sub,nrow = 2) +  
    #plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom",
                                                                                legend.text = element_text(size = 20),
                                                                                legend.title = element_text(size = 20),
                                                                                plot.tag = element_text(size = 30))

figure_2

ggsave("nutr_sub.png", figure_2,
       width = 20, height = 12, dpi = 600)

#B weighted over all imputations----
#0.set-up----
#load data
load(paste0(path, "data_weight_allvisits.Rdata"))

#1. datacleaning----
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

#filter these ID's (defined in part A) from dat_weight_all visits
poor_npcr <- data_weight_allvisits %>%
    filter(id %in% npcr_ids)

poor_bmi <- data_weight_allvisits %>%
    filter(id %in% bmi_ids)

poor_lti <- data_weight_allvisits %>%
    filter(id %in% lti_ids)

poor_sci <- data_weight_allvisits %>%
    filter(id %in% sci_ids)

#make lists for each imputation
imp_list_npcr <- poor_npcr %>%
    group_split(.imp)

imp_list_bmi <- poor_bmi %>%
    group_split(.imp)

imp_list_sci <- poor_sci %>%
    group_split(.imp)

imp_list_lti<- poor_lti %>%
    group_split(.imp)

#2. LMM npcr----
fits_poor_npcr <- map(imp_list_npcr, fit_wemix, outcome = "npcr")

result_poor_npcr <- tidy_we_mix_results(fits_poor_npcr)
#pool
pool_poor_npcr <- rubins_rules_lmm(result_poor_npcr, m = 10)
#remove visit 45 as there is only 1 group left at that point
pool_poor_npcr <- pool_poor_npcr %>%
    filter(term != "visit_time45")

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_npcr <- pool_poor_npcr %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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
    #"visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_npcr %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_npcr)) %>%
    pull(m)

mean_baseline_hdf <- poor_npcr %>%
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
theta <- pool_poor_npcr$theta
ll <- pool_poor_npcr$ll
ul <- pool_poor_npcr$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_npcr$term
names(ll) <- pool_poor_npcr$term
names(ul) <- pool_poor_npcr$term
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

plots_sub_sens <- list()
#plot the predicted values            
plots_sub_sens[[1]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "nPCR (g/kg/day)", y_limits = c(0.3, 1.6))
plots_sub_sens[[1]]


#3.0 SUBGROUP LMM BMI----
fits_poor_bmi <- map(imp_list_bmi, fit_wemix, outcome = "bmi")

result_poor_bmi <- tidy_we_mix_results(fits_poor_bmi)
#pool
pool_poor_bmi <- rubins_rules_lmm(result_poor_bmi, m = 10)

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_bmi <- pool_poor_bmi %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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
    #"visit_time45", "visit_time45:group1",
    "visit_time6", "visit_time6:group1",
    "visit_time9", "visit_time9:group1"
)

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_bmi %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_bmi)) %>%
    pull(m)

mean_baseline_hdf <- poor_bmi %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_bmi)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline bmi value to the corresponding mean value of the group
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

#next, we make vectors from theta, lower and upper limit from pool we mix
theta <- pool_poor_bmi$theta
ll <- pool_poor_bmi$ll
ul <- pool_poor_bmi$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_bmi$term
names(ll) <- pool_poor_bmi$term
names(ul) <- pool_poor_bmi$term
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

#lower then 10 not possible
plot_data_wemix$ll <- if_else(plot_data_wemix$ll <10, 10, plot_data_wemix$ll)

expr_bmi <- expression(BMI (kg/m^2))
#plot the predicted values            
plots_sub_sens[[2]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = expr_bmi, y_limits = c(10, 40))
plots_sub_sens[[2]]




#4.0 SUBGROUP LMM LTI----
fits_poor_lti <- map(imp_list_lti, fit_wemix, outcome = "lti")

result_poor_lti <- tidy_we_mix_results(fits_poor_lti)
#pool
pool_poor_lti <- rubins_rules_lmm(result_poor_lti, m = 10)

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_lti <- pool_poor_lti %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_lti %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_lti)) %>%
    pull(m)

mean_baseline_hdf <- poor_lti %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_lti)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline lti value to the corresponding mean value of the group
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
theta <- pool_poor_lti$theta
ll <- pool_poor_lti$ll
ul <- pool_poor_lti$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_lti$term
names(ll) <- pool_poor_lti$term
names(ul) <- pool_poor_lti$term
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
plots_sub_sens[[3]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = expr_lti, y_limits = c(10, 22))
plots_sub_sens[[3]]


#5.0 SUBGROUP LMM SCI----
fits_poor_sci <- map(imp_list_sci, fit_wemix, outcome = "sci")

result_poor_sci <- tidy_we_mix_results(fits_poor_sci)
#pool
pool_poor_sci <- rubins_rules_lmm(result_poor_sci, m = 10)

selected_terms <- paste0("visit_time", seq(3, 45, by = 3))

table_poor_sci <- pool_poor_sci %>%
    filter(term %in% selected_terms) %>%
    select(term, theta, ll, ul) %>%
    arrange(as.numeric(str_extract(term, "\\d+"))) %>%
    mutate(
        across(
            c(theta, ll, ul),
            ~ round(.x, 3)
        )
    )

#create basic matrix with colnames in the exact same order as the results from pool_we_mix terms
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

#we want the results split for each time point and each group to be able to plot it
#so our dataframe should start with this
plot_data_wemix <- expand.grid(
    visit_time =(c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)),
    group = c(0,1))
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
mean_baseline_hd  <- poor_sci %>%
    filter(group == 0) %>%
    summarise(m = mean(baseline_sci)) %>%
    pull(m)

mean_baseline_hdf <- poor_sci %>%
    filter(group == 1) %>%
    summarise(m = mean(baseline_sci)) %>%
    pull(m)

#residual and id.intercept are variance terms, not coefficients, they are not relevant for prediction so we set them to 0
#the intercept applies for everyone so that should be set to 1
#furthermore, we set the baseline sci value to the corresponding mean value of the group
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
theta <- pool_poor_sci$theta
ll <- pool_poor_sci$ll
ul <- pool_poor_sci$ul
#we match these numbers with the corresponding names of the coefficeints
names(theta) <- pool_poor_sci$term
names(ll) <- pool_poor_sci$term
names(ul) <- pool_poor_sci$term
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
plots_sub_sens[[4]] <- plot_lmm_nutrition_weight(plot_data_wemix, y_lab = "SCI (mg/kg/day)", y_limits = c(10, 25))
plots_sub_sens[[4]]

#COMBINE SUBGROUPS in 1 plot----
figure_sens <- wrap_plots(plots_sub_sens,nrow = 2) +  
    #plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect") & plot_annotation(tag_levels = "A")& theme(legend.position = "bottom", plot.tag=element_text(size=30))

figure_sens

ggsave("nutr_sub_sens.png", figure_sens,
       width = 20, height = 12, dpi = 600)


#C unweighted----
#0. set-up----
load(paste0(path, "data_unweighted.Rdata"))
#1.datacleaning----
data_unweighted <- set_visit_time(data_unweighted, as_factor = TRUE)
#for all analyses, status can be 2 or 0 there is no seperate competing risk in our analyses 
#define npcr scale, so you can perform linear cox for increase per 0.1 in npcr
#we want the hazard for npcr for each 0.1 increase in npcr, so we multiply exposure by 10 by creating npcrscale
data_unweighted <- set_endpt_and_npcrscale(data_unweighted)

#define baseline also for CRP (and Kt/V for other publication)
data_unweighted <- data_unweighted %>%
    group_by(.imp, id)%>%
    mutate(
        baseline_crp_pre_mgdl = crp_pre_mgdl[visit==0],
        baseline_ktv_pre = ktv_pre[visit == 0]
    )%>%
    ungroup()

#filter these ID's (defined in part A) from dat_weight_all visits
poor_npcr <- data_unweighted%>%
    filter(id %in% npcr_ids)

poor_bmi <- data_unweighted %>%
    filter(id %in% bmi_ids)

poor_lti <- data_unweighted %>%
    filter(id %in% lti_ids)

poor_sci <- data_unweighted %>%
    filter(id %in% sci_ids)

#make lists for each imputation
imp_list_npcr <- poor_npcr %>%
    group_split(.imp)

imp_list_bmi <- poor_bmi %>%
    group_split(.imp)

imp_list_sci <- poor_sci %>%
    group_split(.imp)

imp_list_lti<- poor_lti %>%
    group_split(.imp)

#2. LMM npcr----
#apply LMM function each set in the imputation list
fits_npcr_unweight <- map(imp_list_npcr, fit_mixed_model, outcome = "npcr")
#pool them using rubins rules (from package mice)
pooled_npcr_unweight <- pool(fits_npcr_unweight)

pred_list_unweight <- map(fits_npcr_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
#turn into dataframe
preds_unweight <- bind_rows(pred_list_unweight, .id = ".imp")
#define length of imputations
m <- length(fits_npcr_unweight)                        
#manually pool predictions
pooled_pred_unweight <- pooled_pred(preds_unweight) 

plot_unweight <- list()
plot_unweight[[1]] <- plot_lmm_nutrition(pooled_pred_unweight, y_lab = "nPCR (g/kg/day)", y_limits = c(0.5, 1.8))
plot_unweight[[1]]

#3. LMM BMI----
fits_bmi_unweight <- map(imp_list_bmi, fit_mixed_model, outcome = "bmi")
#pool across imputations
pooled_bmi_unweight <- pool(fits_bmi_unweight)
#predict
pred_list_bmi_unweight <- map(fits_bmi_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
#make to dataframe
preds_bmi_unweight <- bind_rows(pred_list_bmi_unweight, .id = ".imp") 
#iterations
m <- length(fits_bmi_unweight)                            
#manually pool predicted values. 
pooled_pred_bmi_unweight <- pooled_pred(preds_bmi_unweight)
expr_bmi <- expression(BMI (kg/m^2))
plot_unweight[[2]] <- plot_lmm_nutrition(pooled_pred_bmi_unweight, y_lab = expr_bmi, y_lim = c(10, 30))
plot_unweight[[2]]

#4.0 LTI----
fits_lti_unweight <- map(imp_list_lti, fit_mixed_model, outcome = "lti")
#pool 
pooled_lti_unweight <- pool(fits_lti_unweight)

#predict
pred_list_lti_unweight <- map(fits_lti_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
preds_lti_unweight <- bind_rows(pred_list_lti_unweight, .id = ".imp")  
m <- length(fits_lti_unweight)                            

pooled_pred_lti_unweight <-pooled_pred(preds_lti_unweight) 
expr_lti <- expression (LTI (mg/m^2))
plot_unweight[[3]] <- plot_lmm_nutrition(pooled_pred_lti_unweight, y_lab = expr_lti, y_lim = c(10, 22))
plot_unweight[[3]]

#5.0 sci----
fits_sci_unweight <- map(imp_list_sci, fit_mixed_model, outcome = "sci")
#pool
pooled_sci_unweight <- pool(fits_sci_unweight)
#predict
pred_list_sci_unweight <- map(fits_sci_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
preds_sci_unweight <- bind_rows(pred_list_sci_unweight, .id = ".imp")  
m <- length(fits_sci_unweight)                             

pooled_pred_sci_unweight <- pooled_pred(preds_sci_unweight)

plot_unweight[[4]] <- plot_lmm_nutrition(pooled_pred_sci_unweight, y_lab = "SCI (mg/kg/day)", y_lim = c(15, 25))
plot_unweight[[4]]


figure_sens_unweight <- wrap_plots(plot_unweight,nrow = 2) +  
    #plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect") & plot_annotation(tag_levels = "A")& theme(legend.position = "bottom", plot.tag = element_text(size = 30))

figure_sens_unweight

ggsave("nutr_sub_unweight.png", figure_sens_unweight,
       width = 20, height = 12, dpi = 1200)
