#imputation for CONVINCE longitudinal data, for both nutritional and UFR analyses
#by Sanne Roos
#last updated on 5 dec 2025


#0. set-up ----
##load packages----
pacman::p_load ("rio",     #to recognize other files and load them in R
                "conflicted", #if a function is in more packages it makes you choose which one to use
                "tidyverse", #to work with nice syntaxes
                "mice", #for imputation and missing data tests
                "VIM", #visualization and imputation of missing values
                "miceadds", #for different imputation method for longitudinal data
                "readxl" #read excel
)

###resolve package conflicts----
conflicts_prefer(dplyr::filter) #between dplyr and stats
conflicts_prefer(dplyr::lag) #between dplyr and stats

####set working directory
path <- "/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/"

#####read database----
before_imp_data <- read_rds(file ="/Users/sroos6/Documents/Projecten/HDF pooling trials/data/convince_nutrition_1.rds")

#With comorbidity variables
comorb <- read_excel("/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/data/Archief 2/CONVINCE_export_20230417.xlsx")
#Define baseline comorbidity variables of interest in a vector
comorb_var <- c("Participant Id", "SCR_MH_AP_OCCUR",  "SCR_MH_MYO_OCCUR",            "SCR_MH_PTCA_OCCUR",        "SCR_MH_CABG_OCCUR",            
                 "SCR_MH_PACE_OCCUR",           "SCR_MH_DEFIB_OCCUR",         "SCR_MH_CHF_OCCUR",         "SCR_MH_AFIB_OCCUR",          "SCR_MH_TIA_OCCUR",            
                 "SCR_MH_CVA_OCCUR",              "SCR_MH_CEA_OCCUR",              "SCR_MH_INTCLAU_OCCUR",          "SCR_MH_PTA_OCCUR",              "SCR_MH_BYPLL_OCCUR",           
                 "SCR_MH_LLIMBAMP_OCCUR",         "SCR_MH_AAA_OCCUR" ,             "SCR_MH_DOTSTEN_OCCUR",          "SCR_MH_DIAB_OCCUR" ,            "SCR_MH_CAN_OCCUR" ,            
                 "SCR_MH_COPD_OCCUR" )

#add these to our main database. Note that in in the comorbidity database, id = Participant id
before_imp_data_merged <- before_imp_data %>%
    left_join(
        #take from comorbidity the participant id and all the variables of interest as defined in comorb_var
        comorb %>% select(`Participant Id`, all_of(comorb_var)), 
        by = c("id" = "Participant Id") # match on id
    )
#give missing visit dates a flag to remove them after imputation
#if you do not do this longitdunal imputation fails for ID's with only one observation
before_imp_data_merged <- before_imp_data_merged %>%
    filter(!is.na(visit_date))

#skip visit 16 as there is only one patient there
before_imp_data_merged <- before_imp_data_merged %>%
    filter(visit != 16)

#select variables of interest for imputation process
before_imp_data_merged <- before_imp_data_merged %>%
    select("id","name","country","group","sex","age","visit","visit_date","dial_vintage",
          "height","sbp_pre","dbp_pre","uf_vol","dial_time", "sbp_post", "dbp_post",
          "weight_post", "bsa_dub","ktv_pre","blood_flow","hb_pre_mmoll",
          "creat_pre_umoll","urea_pre_mgdl","crp_pre_mgdl", "t_start", "t_stop", "endpt",
          "pth_pre_pmoll","na_pre_mmoll","k_pre_mmoll","ca_pre_mmoll","phos_pre_mmoll",
          "mg_pre_mmoll","residual_urine_out", "hdf_convol", "end_eot_reason_recode_full", "ran_date", "end_eos_reason_recode_full",
          SCR_MH_AP_OCCUR:SCR_MH_COPD_OCCUR)

save(before_imp_data_merged, file = paste0(path, "before_imp_data_merged.Rdata"))


#1.0 additional datacleaning----
#dialysis vintage does not change in time so we fill it for each row. (we do not want imputation there)
before_imp_data_merged <- before_imp_data_merged %>%
    group_by(id) %>%
    fill(dial_vintage, .direction = "downup") %>%
    mutate(
        #set group to 1 and 0, which is more easy for imputation model
        group = case_when(
            group == "High-dose hemodiafiltration" ~ 1,
            group == "High-flux hemodialysis" ~ 0), 
        group = as.factor(group), 
        name = str_replace_all(name, "\\s+", ""),#remove spaces from name
        end_eot_reason_recode_full = as.factor(end_eot_reason_recode_full),
        end_eos_reason_recode_full = as.factor(end_eos_reason_recode_full)
    )%>%
    ungroup()
#we set hdf convective volume to 0 for everyone on HD because this does not apply to them and we do not want imputation there
before_imp_data_merged[["hdf_convol"]] <- as.numeric(before_imp_data_merged[["hdf_convol"]])
before_imp_data_merged[["hdf_convol"]] <- if_else(before_imp_data_merged[["group"]] == 0, 0, before_imp_data_merged[["hdf_convol"]])

#2 imputation process----
#2.0 get prediction matrix
mat_prd <- mice(before_imp_data_merged, maxit = 0)[["predictorMatrix"]]

#2.1 Define variables that do not have to be imputed
vec_nimp <- c(# Auxiliary variables
    "id",  "country", "group", "visit", #identifier or time indicator variables
    # Dates
    "visit_date", "ran_date", "t_start", "t_stop", "endpt",
    #other predictive variables
    "sex" ,                 "age",                   "dial_vintage",         "height", 
    #comorbidity
    "SCR_MH_AP_OCCUR",       "SCR_MH_MYO_OCCUR",     
    "SCR_MH_PTCA_OCCUR",    "SCR_MH_CABG_OCCUR",     "SCR_MH_PACE_OCCUR", 
    "SCR_MH_DEFIB_OCCUR",    "SCR_MH_CHF_OCCUR",      "SCR_MH_AFIB_OCCUR",     "SCR_MH_TIA_OCCUR",     
    "SCR_MH_CVA_OCCUR",      "SCR_MH_CEA_OCCUR",      "SCR_MH_INTCLAU_OCCUR",  "SCR_MH_PTA_OCCUR",  
    "SCR_MH_BYPLL_OCCUR",    "SCR_MH_LLIMBAMP_OCCUR", "SCR_MH_AAA_OCCUR",     
    "SCR_MH_DOTSTEN_OCCUR",  "SCR_MH_DIAB_OCCUR",     "SCR_MH_CAN_OCCUR",      "SCR_MH_COPD_OCCUR",
    #longitudinal variables
    "sbp_pre",    "dbp_pre",                                       "bsa_dub",              
     "blood_flow",                       "hb_pre_mmoll",            "na_pre_mmoll" ,  
    "k_pre_mmoll",                    "ca_pre_mmoll",           "phos_pre_mmoll", 
    "mg_pre_mmoll",           "crp_pre_mgdl",        "pth_pre_pmoll",             
     "sbp_post", "dbp_post", "residual_urine_out", "hdf_convol")

#define the dates
date_vars <- c(
   "visit_date", "ran_date")
#define text variables
text_vars <- c("end_eot_reason_recode_full", "end_eos_reason_recode_full")

# Define variables for longitudinal imputation, we also need uf_vol and dial_time to calculate ufr
#we need urea pre mgdl to recalculate BUN, and after that we need ktv pre to calculate npcr
#for UFR we need dry weight, uf volume and dial time
vec_limp <- c(
     "urea_pre_mgdl", "ktv_pre", "creat_pre_umoll", "weight_post", "height", "uf_vol", "dial_time")

# Adjust predictor matrix accordingly
mat_prd[vec_nimp, ] <- 0

#make sure dates and text variables are not used to impute
mat_prd[, colnames(mat_prd) %in% date_vars | colnames(mat_prd) %in% text_vars] <- 0

# Set the cluster variable for longitudinal imputation
mat_prd[, "id"] <- -2

# Get methods vector
vec_mtd <- mice(before_imp_data_merged, maxit = 0)[["method"]]

# Change variables to be longitudinally imputed to longitudinal imputation
vec_mtd <- if_else(names(vec_mtd) %in% vec_limp, "2l.pmm", vec_mtd)

#use no method for the dates
vec_mtd[names(vec_mtd) %in% date_vars| names(vec_mtd) %in% text_vars] <- ""

# Reset names of methods vector
names(vec_mtd) <- colnames(mat_prd)

# Start imputation, 10 datasets(m), 50 iterations (maxit)
lst_imp <- mice(before_imp_data_merged,
                m = 10,
                maxit = 50,
                method = vec_mtd,
                predictorMatrix = mat_prd,
                seed = 1) # set seed to make sure it can be repeated

#NOTE: running this imputation results in several warnings regarding scaling, this does not effect the quality of imputation result, we can ignore this
#We also get message that model fails to converge in some of the iterations. 
#Plotting the result shows random plots for imputation, means we are satisfied with the result:
plot(lst_imp)
# Save imputation object
save(lst_imp, file = paste0(path, "imputation_object.Rdata"))

#3.0 Finalise imputed data ----
#Load imputation object
load(paste0(path, "imputation_object.Rdata"))

#convert to long format dataframe using complete
dat_imputed <- complete(lst_imp, 
                        action = "long",
                        include = TRUE) %>%
    group_by(id) %>%
    # Calculate extra variables
    mutate(# BMI
        bmi = weight_post / (height / 100) ^ 2,
        bun = urea_pre_mgdl / 2.1428,
        npcr = (bun / (25.8 + (1.15 * ktv_pre) + 56.4/ktv_pre)) + 0.168,
        sci = (16.21 + if_else(sex == 0, 1.12, 0) - (0.06 *age) - (0.08 * as.numeric(ktv_pre)) + (0.009 * as.numeric(creat_pre_umoll))),
        lti = ((sci * weight_post*0.029)+7.38)/((height/100)^2),
        ultra_fil_rate = uf_vol / (weight_post*dial_time/60)
    ) %>%
    ungroup()
#save
save(dat_imputed, file = paste0(path, "imputed_data_convince.Rdata"))
