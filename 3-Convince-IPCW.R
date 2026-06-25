#IPCW
#by Sanne Roos
#last updated on march 19 2026

#0. set-up ----
#load packages
pacman::p_load ("rio",     #to recognize other files and load them in R
                "conflicted", #if a function is in more packages it makes you choose which one to use
                "tidyverse", #tidy data
                "ggplot2" #make plots
)

##resolve package conflicts
conflicts_prefer(dplyr::filter) #between dplyr and stats
conflicts_prefer(dplyr::lag) #between dplyr and stats
conflicts_prefer(magrittr::set_names) # between magrittr and purrr
###set working directory
path <- "/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/"
####load data
load(paste0(path, "imputed_data_convince.Rdata"))

#source functions
walk(list.files(paste0(path, "R/funs")), ~ source(paste0(path, "R/funs/", .x)))

#1.0 prepare for censoring----
#run on imputed data, skip the unimputed first dataset
dat_imputed <- dat_imputed %>%
    filter(.imp !=0)

#the primary analysis for CONVINCE was intention to treat, meaning they kept following the patient even after NTx
#make sure censor is 1 at the last row of each individual in each imputed dataset, 
#so censor if they stopped before study end due to informative censoring reason
#informative censoring reasons are:
#stopping dialysis, kidney transplant, changed modality, patient decision, other or dying.
#all these lead to stopping of nutritional measures (or ultrafiltration measures (different project) and it is likely their nutritional status changes due to these events
data_unweighted <- dat_imputed %>%
    group_by(.imp, id) %>%
    mutate(
      inf_cens = case_when(
          endpt == 1 | endpt == 2 ~ 1, #both death (1) or competing risk (coded as 2) are informative
          endpt == 0 ~ 0), #uninformative censoring
        #set baseline values to later use in the linear mixed model
        baseline_npcr = npcr[visit == 0],
        baseline_ultra_fil_rate = ultra_fil_rate[visit==0],
        baseline_bmi = bmi[visit ==0],
        baseline_sci = sci[visit ==0],
        baseline_lti = lti[visit == 0],
        baseline_weight = weight_post[visit == 0]
    ) %>%
    ungroup()

#store in between as data without weighting
save(data_unweighted, file = paste0(path, "data_unweighted.Rdata"))

#1.1 check positivity for each visit in each imputation----
#it seems that for visit 15 in most imputations there is no censoring. 
#Futhermore, it seems that at each visit, there may be too little events to make a meaningful censoring model for each visit.
#we continue and plot the result later in this script
by(data_unweighted$endpt, list(data_unweighted[["visit"]], data_unweighted[[".imp"]]), table)


#2.0 ipcw combining all visits ----
#make a vector to loop over each imputation
imps <- sort(unique(data_unweighted[[".imp"]]))

#with map you can apply a function to each element of a vector and put it in a list,
#in this case, the the generalized linear censoring model for each imputation
fit_list <- map(imps, glm_for_weighting_across_imp)
names(fit_list) <- paste0("imp_", imps) #to inspect fit list we put imp before it

#give predictions and put in a list
ps_list <- map(imps, predict_weight)
names(ps_list) <- paste0("imp_", imps)

#with map2 you can iterate over two arguments at the time, in this case imps and ps_list
#ps_list is a list with probabilities for each imputation
#so essentially, what we do here:
#map2(imps, ps_list, function(i, ps) { ... }) returns dat_imp for ach combination of i and ps:
#Iteration 1: i = 1, ps = ps_list[[1]]
#dat_imp = all rows of .imp == 1
#Add columns (ps_ipcw, ipcw)
#Return dat_imp
#Iteration 2: i = 2, ps = ps_list[[2]] etcetera..
data_weight_allvisits <- bind_rows( #turn into dataframe
    map2(imps, ps_list, prob_to_ipcw))

#we cap the weights because most weights around 1 but maximum is 190. 
# we cap at 99 percentile, because we have few weights of 800. 99th percentile is at 34.9
cap_99 <- quantile(data_weight_allvisits[["ipcw"]], 0.99)
#if ipcw is larger than cap99, take cap 99.
data_weight_allvisits <- data_weight_allvisits %>%
    mutate(ipcw = pmin(ipcw, cap_99))

#we try to make a density plot to see if weighing went well.
unweightplot_weight_all_visits <- unweight_plot_ipcw(data_weight_allvisits)
unweightplot_weight_all_visits

weightplot_weight_all_visits <- weight_plot_ipcw(data_weight_allvisits)
weightplot_weight_all_visits

#we can also plot this for each visit on example imputation 1
#example on first imputation for plotting
example <- data_weight_allvisits %>%
    filter(.imp == 1)
plots <- list()
plots[[1]] <- plot_ipcw_visit_unweighted(data = example) + theme(legend.position = "none")
plots[[1]]

plots[[2]] <- plot_ipcw_visit_weighted(data = example) + theme(legend.position = "none")
plots[[2]]

figure_1 <- wrap_plots(plots, ncol = 1) +  
    plot_annotation(tag_levels = "A") +
    theme(legend.position = "bottom")
figure_1
#save it
ggsave("plot_weighted&unweighted.png", figure_1, width= 8, height = 10, dpi = 1200)

#save inbetween, the weights for each imputation, not fitted per visit
save(data_weight_allvisits, file = paste0(path, "data_weight_allvisits.Rdata"))



#3.0 now apply for each imputation and visit seperately----
visits<- 0:13
fit_list3 <- map(imps, glm_for_weighting_total)
#set the first names to imp_
names(fit_list3) <- paste0("imp_", imps)
#set names within the fit_list 3 to the right visits
fit_list3 <- map(fit_list3, function(.x) set_names(.x, visits))

ps_list <- map(imps, predict_weight_total)
#now we have a list with 10 dataframes, we want to turn that into one single dataframe
ps_list_combined <- bind_rows(ps_list)
#now we calculat ipcw from he ps
data_weight_allvisits_imp <- ps_list_combined %>%
    mutate(
        ps_ipcw = ps,
        ipcw = case_when(
            inf_cens == 1 ~ 1 / (1 - ps_ipcw),  
            inf_cens == 0 ~ 1 / ps_ipcw 
        )
    )

# we cap at 99 percentile, because we have few weights of 430. 99th percentile is at 30.7
cap_99 <- quantile(data_weight_allvisits_imp[["ipcw"]], 0.99)
#if ipcw is larger than cap99, take cap 99.
data_weight_allvisits_imp <- data_weight_allvisits_imp %>%
    mutate(ipcw = pmin(ipcw, cap_99))

#combine with visits 14 and 15 again
#left_visits <- data_unweighted %>%
    #filter(visit == 14 | visit == 15) %>%
    #mutate(
        #ps = 1,
        #ps_ipcw = 1,
        #ipcw = 1)

#data_weight_allvisits_imp <- rbind(data_weight_allvisits_imp, left_visits) %>%
    #arrange(id, .imp, visit)


save(data_weight_allvisits_imp, file = paste0(path, "data_weight_allvisits_imp.Rdata"))

border_theme <- theme(
    plot.background = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
    ),
    plot.margin = margin(5, 5, 5, 5)
)
plots <- list()
plots[[1]] <- plot_ipcw_visit_unweighted(data = data_weight_allvisits_imp) + border_theme + theme(legend.position = "none")
plots[[1]]

plots[[2]] <- plot_ipcw_visit_weighted(data = data_weight_allvisits_imp) + border_theme+ theme(legend.position = "none")
plots[[2]]

figure_1 <- wrap_plots(plots, ncol = 2) +  
    plot_annotation(tag_levels = "A") +
    theme(legend.position = "bottom")
figure_1

ggsave("plot_weighted&unweighted.png", figure_1, width= 15, height = 15, dpi = 1200)
