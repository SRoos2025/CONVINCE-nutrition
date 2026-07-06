#nutrition analyses LMM and Cox
#script by S. Roos 
#last updated 18-06-2026

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
               "lme4", #linear mixed models
               "ggeffects", # to predict values from linear mixed models
               "broom.mixed", # supports mixed models, returns tidy output
               "ggplot2",# to plot the predicted values) # Printing tables to word ) 
               "survival",#to use cox and penalized spline) 
               "gtsummary",#to make a formatted table )
               "splines")#for robust standard error with weighted linear mixed model
#solve conflicts 
conflicts_prefer(dplyr::filter) # Between dplyr & stats 
conflicts_prefer(dplyr::select) # Between dplyr & MASS
conflicts_prefer(lmerTest::lmer) #between lmerTest and lme4

#set working directory 
path <- "/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/"

setwd("/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/")

#load data. 
#we perform 3 analyses, unweighted, weighted across all imputations and weighted across all imputations AND visits
#in this script, only the unweighted
#unweighted
load(paste0(path, "data_unweighted.Rdata"))

# Source relevant functions
walk(list.files(paste0(path, "R/funs")), ~ source(paste0(path, "R/funs/", .x)))

#1.0 datacleaning----
#make visit a factor and #set fixed times for each visit, take the mean of each visit
#you can either use the factor or not in the function, if you want time continuous for the estimates, set to FALSE
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

#make a list for each dataframe of each imputed dataset
imp_list_unweight <- data_unweighted %>% 
    group_split(.imp)

#1.1LMM unweighted----
##1.1.2 LMM nPCR unweighted----
#apply LMM function each set in the imputation list
fits_npcr_unweight <- map(imp_list_unweight, fit_mixed_model, outcome = "npcr")
#pool them using rubins rules (from package mice)
pooled_npcr_unweight <- pool(fits_npcr_unweight)
#we want a summary including confidence intervals
summary <- summary(pooled_npcr_unweight, conf.int = TRUE)
#turn into a table
knitr::kable(
    summary,
    digits = 4,
    caption = "Pooled linear mixed-effects model results for NPCR"
)

##1.1.3 unweighted LMM plot of npcr predicted values----
#make a list so other nutritional variables can also be plotted in there
plots <- list()
#create predicted values per group
pred_list_unweight <- map(fits_npcr_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
#turn into dataframe
preds_unweight <- bind_rows(pred_list_unweight, .id = ".imp")
#define length of imputations
m <- length(fits_npcr_unweight)                        
#manually pool predictions
pooled_pred_unweight <- pooled_pred(preds_unweight) 

plots[[1]] <- plot_lmm_nutrition(pooled_pred_unweight, y_lab = "nPCR (g/kg/day)", y_limits = c(0.5, 1.8))
plots[[1]]
###1.1.4 LMM for BMI----
#apply lmm
fits_bmi_unweight <- map(imp_list_unweight, fit_mixed_model, outcome = "bmi")
#pool across imputations
pooled_bmi_unweight <- pool(fits_bmi_unweight)
#we want a summary including confidence intervals
summary_bmi_unweight<- summary(pooled_bmi_unweight, conf.int = TRUE)
#show table
knitr::kable(
    summary_bmi_unweight,
    digits = 3,
    caption = "Pooled linear mixed-effects model results for BMI"
)
###1.1.5 plot predicted values BMI----
#created predicted values per group
pred_list_bmi_unweight <- map(fits_bmi_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
#make to dataframe
preds_bmi_unweight <- bind_rows(pred_list_bmi_unweight, .id = ".imp") 
#iterations
m <- length(fits_bmi_unweight)                            
#manually pool predicted values. 
pooled_pred_bmi_unweight <- pooled_pred(preds_bmi_unweight)
expr_bmi <- expression(BMI (kg/m^2))
plots[[2]] <- plot_lmm_nutrition(pooled_pred_bmi_unweight, y_lab = expr_bmi, y_lim = c(20, 40))
plots[[2]]
    
####1.1.6 LMM for LTI----
#apply LMM function on each set in the imputation list
fits_lti_unweight <- map(imp_list_unweight, fit_mixed_model, outcome = "lti")
#pool 
pooled_lti_unweight <- pool(fits_lti_unweight)
#we want a summary including confidence intervals
summary_lti_unweight <- summary(pooled_lti_unweight, conf.int = TRUE)
#summarize in table
knitr::kable(
    summary_lti_unweight,
    digits = 3,
    caption = "Pooled linear mixed-effects model results for LTI"
)

####1.1.76 plot predicted values LTI----
pred_list_lti_unweight <- map(fits_lti_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
preds_lti_unweight <- bind_rows(pred_list_lti_unweight, .id = ".imp")  
m <- length(fits_lti_unweight)                            

pooled_pred_lti_unweight <-pooled_pred(preds_lti_unweight) 
expr_lti <- expression (LTI (mg/m^2))
plots[[3]] <- plot_lmm_nutrition(pooled_pred_lti_unweight, y_lab = expr_lti, y_lim = c(15, 22))
plots[[3]]
    

#####1.1.8 LMM for SCI----
#apply LMM function on each set in the imputation list
fits_sci_unweight <- map(imp_list_unweight, fit_mixed_model, outcome = "sci")
#pool
pooled_sci_unweight <- pool(fits_sci_unweight)
#we want a summary including confidence intervals
summary_sci_unweight <- summary(pooled_sci_unweight, conf.int = TRUE)
#put in table
knitr::kable(
    summary_sci_unweight,
    digits = 3,
    caption = "Pooled linear mixed-effects model results for SCI"
)

#####1.1.9 predict and plot from LMM SCI----
pred_list_sci_unweight <- map(fits_sci_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
preds_sci_unweight <- bind_rows(pred_list_sci_unweight, .id = ".imp")  
m <- length(fits_sci_unweight)                             

pooled_pred_sci_unweight <- pooled_pred(preds_sci_unweight)

plots[[4]] <- plot_lmm_nutrition(pooled_pred_sci_unweight, y_lab = "SCI (mg/kg/day)", y_lim = c(15, 25))
    plots[[4]]

    #1.1.1.1 combine LMM plots in one----
    figure_unweight <- wrap_plots(plots,nrow = 2) +  
        #plot_annotation(tag_levels = "A") +
        plot_layout(guides = "collect") & plot_annotation(tag_levels = "A")& theme(legend.position = "bottom", plot.tag = element_text(size = 30))
    
    figure_unweight
    
    ggsave("nutr_unweight.png", figure_unweight,
           width = 20, height = 12, dpi = 1200)

# marge between plos
figure_1 <- figure_1 & theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2), legend.position = "bottom",  legend.key.size = unit(0.5, "cm"),  legend.text = element_text(size = 14))  
figure_1
#save
ggsave("lmm_unweighted_nutr.png", figure_1,
       width = 14, height = 10, dpi = 1200)

#LMM CRP----
fits_crp_unweight <- map(imp_list_unweight, fit_mixed_model, outcome = "crp_pre_mgdl")
#pool
pooled_crp_unweight <- pool(fits_crp_unweight)
#we want a summary including confidence intervals
summary_crp_unweight <- summary(pooled_crp_unweight, conf.int = TRUE)
#plot predicted CRP
pred_list_crp_unweight <- map(fits_crp_unweight, ~ ggpredict(.x, terms = c("visit_time", "group")) %>% as.data.frame())
preds_crp_unweight <- bind_rows(pred_list_crp_unweight, .id = ".imp")  
m <- length(fits_crp_unweight)                           

pooled_pred_crp_unweight <-pooled_pred(preds_crp_unweight) 

plot_crp <- plot_lmm_nutrition(pooled_pred_crp_unweight, y_lab = "CRP(mg/dL)", y_lim = c(0, 5))
plot_crp

ggsave("crp_unweighted.png", plot_crp,
       width = 14, height = 10, dpi = 1200)

#3.0 linear cox unweighted---- 
#for association nutrition and mortality
##3.0 prepare data----


#make again implist with now also variable npcr_scale
imp_list_unweight <- data_unweighted %>% 
    group_split(.imp)

###3.1 nPCR----
#run cox function on each imputation
fits_npcr_cox_unweight <- map(imp_list_unweight, \(x)cox_nutr(x, "npcr_scale", robust = FALSE, weights = FALSE))
#pool each imputation
pooled_npcr_cox_weight <- pool(fits_npcr_cox)
#put in in a tidy table
result_npcr <- broom::tidy(pooled_npcr_cox, exponentiate = TRUE, conf.int = TRUE)

####3.2 BMI----
#\(x) anonymous function, we do not re-use it but immediately define it
fits_bmi_cox <- map(imp_list, \(x)cox_nutr(x, "bmi"))
#pool
pooled_bmi_cox <- pool(fits_bmi_cox)
#we want hazard ratio so exponentiate = TRUE converts this from estimate to exponent
result_bmi <- broom::tidy(pooled_bmi_cox, exponentiate = TRUE, conf.int = TRUE)

#####3.3 SCI----
#fit cox on each impuation
fits_sci_cox <- map(imp_list, \(x)cox_nutr(x, "sci"))
#pool
pooled_sci_cox <- pool(fits_sci_cox)
#we want hazard ratio so exponentiate = TRUE converts this from estimate to exponent
result_sci <- broom::tidy(pooled_sci_cox, exponentiate = TRUE, conf.int = TRUE)

######3.4 LTI----
#fit cox on each imputation
fits_lti_cox <- map(imp_list, \(x)cox_nutr(x, "lti"))
#pool
pooled_lti_cox <- pool(fits_lti_cox)
#we want hazard ratio so exponentiate = TRUE converts this from estimate to exponent
result_lti <- broom::tidy(pooled_lti_cox, exponentiate = TRUE, conf.int = TRUE)

#3.5 linear cox weighted----
#weighted per visit AND imputation
fits_npcr_cox_weight <- map(imp_list_weight_imp, \(x)cox_nutr(x, "npcr_scale"))
#extract results and put in dataframe
result_cox_weight <- extract_cox_weight(fits_npcr_cox_weight)
#manually pool using rubins rules
#the function for pooling for lmm works the same for the cox
pooled_result_cox_weight <- rubins_rules_lmm(result_cox_weight, m = 10)


#4.0 normalised spline models----
#check for potential non-linear associations
#we use resource: https://cran.r-project.org/web/packages/survival/vignettes/splines.pdf

##4.1 nPCR unweighted----
#we apply a general function to make cox splines: cox_spline_nutr
#in this function, we must define our nutrition variable, 
#but also the lower and upper bound, because these can not be the lower and upper bound of each imputed dataset
#because if we apply rubin's rules, then the estimates will be for different npcr values in each imputed dataset
#because each imputed dataset has different upper and lower bounds
#so we determine general upper and lower bounds of the entire data_unweighted
fit_spline_npcr_cox <- map(imp_list_unweight, \(x)cox_spline_nutr(x, nutr_var = "npcr", 
                                                         lower_bound = min(data_unweighted[["npcr"]]),
                                                         upper_bound = max(data_unweighted[["npcr"]]),
                                                         weights = FALSE)
)
#turn it from list to dataframe
fit_spline_npcr_cox <- bind_rows(fit_spline_npcr_cox)
#apply rubin's rules
fit_spline_npcr_cox <- rubin_rule_cox_spline(fit_spline_npcr_cox, nutr_var = "npcr", center_at = 1.1)
#make a list for the plots
plots_spline_cox <- list()

#plot it
plots_spline_cox[[1]]<- nutr_flex_plot(fit_spline_npcr_cox, nutr_var = "npcr", x_lab = "nPCR (g/kg/day)", center_val = 1.1)

plots_spline_cox[[1]]
###4.2 BMI----
#we center at 25
fit_spline_bmi_cox <- map(imp_list_unweight, \(x)cox_spline_nutr(x, nutr_var = "bmi", 
                                                        lower_bound = min(data_unweighted[["bmi"]]),
                                                        upper_bound = 50))

#turn it from list to dataframe
fit_spline_bmi_cox <- bind_rows(fit_spline_bmi_cox)
#apply rubin's rules
fit_spline_bmi_cox <- rubin_rule_cox_spline(fit_spline_bmi_cox, nutr_var = "bmi", center_at = 23.9)


#plot it
plots_spline_cox[[2]]<- nutr_flex_plot(fit_spline_bmi_cox, nutr_var = "bmi", x_lab = expr_bmi, center_val = 25, y_limits = c(0, 5), break_min = 10, break_max = 50, breaks = 10)

plots_spline_cox[[2]]

####4.3 LTI----
#for estimated LTI, it is less clear what an "ideal" reference value should be
#therefore, we use median of total group which is 17.2 
fit_spline_lti_cox <- map(imp_list_unweight, \(x)cox_spline_nutr(x, nutr_var = "lti", 
                                                        lower_bound = min(data_unweighted[["lti"]]), upper_bound = max(data_unweighted[["lti"]])))

#turn it from list to dataframe
fit_spline_lti_cox <- bind_rows(fit_spline_lti_cox)
#apply rubin's rules
fit_spline_lti_cox <- rubin_rule_cox_spline(fit_spline_lti_cox, nutr_var = "lti", center_at = 17.2)


#plot it
plots_spline_cox[[3]]<- nutr_flex_plot(fit_spline_lti_cox, nutr_var = "lti", x_lab = expr_lti, center_val = 17.2, break_min = 10, break_max = 40, breaks = 10)
plots_spline_cox[[3]]
#####4.4 SCI----
#for estimated SCI, it is less clear what an "ideal" reference value should be
#therefore, we use the median of the total study population
fit_spline_sci_cox <- map(imp_list_unweight, \(x)cox_spline_nutr(x, nutr_var = "sci", 
                                                        lower_bound = min(data_unweighted[["sci"]]), upper_bound = max(data_unweighted[["sci"]])))

#turn it from list to dataframe
fit_spline_sci_cox <- bind_rows(fit_spline_sci_cox)
#apply rubin's rules
fit_spline_sci_cox <- rubin_rule_cox_spline(fit_spline_sci_cox, nutr_var = "sci", center_at = 19.3)
#plot it
plots_spline_cox[[4]]<- nutr_flex_plot(fit_spline_sci_cox, nutr_var = "sci", x_lab = "SCI (mg/kg/day)", center_val = 19.3, break_min = 15, break_max = 30, breaks = 5)

plots_spline_cox[[4]]
######3.5 Overall plot----
#fuse them together using wrap_plots
figure_1 <- wrap_plots(plots_spline_cox) + plot_annotation(tag_levels = "A")
#check
figure_1 
#save
ggsave("flex_spline_nutr_unweighted.png", figure_1, width= 14, height = 10, dpi = 1200)


