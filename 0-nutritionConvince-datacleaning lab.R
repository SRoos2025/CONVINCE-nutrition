#check variables of CONVINCE study
#by Sanne Roos
#last updated on 18 dec 2025

#0. set-up ----
#load packages
pacman::p_load ("rio",     #to recognize other files and load them in R
                "tidyverse" #tidy data
)

##1.0 import data----
setwd ("/Users/sroos6/Documents/Projecten/HDF pooling trials/data/Archief 2")
data_co <-readRDS("/Users/sroos6/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projecten/10-2025 HDF pooling trials/data/CONVINCE-RScripts-and-Datasets_01Sep2024/2Data/3Cleaned-and-Formatted-Datasets/1Main-Dataset/MainDataset_FINAL_20230417.rds") %>% #read.csv2("CONVINCE_export_20230417.csv") %>%
    #make al variable lower case
    set_colnames(tolower(colnames(.)))%>%
    #we want to remove scr or v00 before the visit
    #To do that, we use string replace. Note: this replaces the first match, string replace_all replaces all matches.
    #In string replace, first state where the string can be found, in this case in colnames(.)
    #Second, define the pattern and then by what to replace.
    #Here we we use reggex. ^ means begin of the string, $ means end of the string.
    #we want the scr to be replaced by v00 because it is baseline, or "screen"
    # with ?<= we look back. So from _ we look back
    #we look back at the beginning so ^, then a v and then there are exactly 2 digits. THIS IS NOT REPLACED. It serves only to recognize
    #that is \d->https://regex101.com/ because in R \ is an expression we escape that by again applying \. We want to replace the _ by . 
    #(because there are other _ in the columnames so we specify the first one and replace by a dot)
    set_colnames(str_replace(str_replace(colnames(.), "^scr_", "v00_"), "(?<=^v\\d{2})_", ".")) %>%
    mutate(
        #From colomn v00.dbp_pre until v17.creat_post_mgdl we want to make it all characters.
        across(v00.visit_date:v17.conmeds_yn, as.character))

#2.0 datacleaning----
#make dataset long in stead of wide
#We want to do this for a selection of columns because for sex and age for example wide format is fine.
#with names_sep, we state that we want different names for the columns and split the column name at . 
#Because it is not in a pipe we must use \\ (see explanation above)
#in names_to we use .value so we dont need to specify values_to anmyore
#it automatically indicates that the corresponding component of the column name defines the name of the output column containing the cell values
#in names_transform, we say that in visit we do not want any numbers anymore so we use parse_number, which strips non-numeric components.
data_long <- pivot_longer(data_co,
                          cols = matches("^v\\d{2}\\..+"),
                          names_sep = "\\.", #everything before the . must become visit, for example v00 becomes visit. Everything after remains what it was with .value
                          names_to = c("visit", ".value"), 
                          names_transform = list(visit = ~ parse_number(.x)) 
)

#save in between
saveRDS(data_long, file ="/Users/sroos6/Documents/Projecten/HDF pooling trials/data/convince_data_long.rds")
#open dataset
#data_long <- read_rds("/Users/sroos6/Documents/Projecten/HDF pooling trials/data/convince_data_long.rds")

#We test how many NA's we have before conversion to numeric, to check this remains similar after conversion to numeric
#creatinin
sum(is.na(data_long$creat_pre_umoll)) #14507
sum(is.na(data_long$creat_pre_mgdl)) #14507
sum(is.na(data_long$creat_pre_mgl)) #14507

#urea
sum(is.na(data_long$urea_pre_mmoll)) #12095
sum(is.na(data_long$urea_pre_gl)) #21198
sum(is.na(data_long$urea_pre_mgdl))#12095

#weight
sum(is.na(data_long$height)) #23124
sum(is.na(data_long$weight_post)) #11920

#now convert to numeric
data_long <- data_long %>%
    #we want the sex variable with the other demographics, befor the visits, for this you can use relocate
    relocate(sex, .before = visit) %>%
    mutate(
        across(age:eot_visit_yn, as.numeric))
#this gives a warning:"NAs introduced by coercion"
#we will check if it is in any of these variables
#creat
sum(is.na(data_long$creat_pre_umoll)) #14507
sum(is.na(data_long$creat_pre_mgdl)) #14507
sum(is.na(data_long$creat_pre_mgl)) #14507

#ureum
sum(is.na(data_long$urea_pre_mmoll)) #12095
sum(is.na(data_long$urea_pre_gl)) #21198
sum(is.na(data_long$urea_pre_mgdl))#12095

#weight
sum(is.na(data_long$height)) #23124
sum(is.na(data_long$weight_post)) #11920

#these NA's have not changed compared to the NA's above, before conversion to numeric
#so we continue
#now show summary and distribution of creat and ureum
hist(data_long$creat_pre_umoll)
summary(data_long$creat_pre_umoll) #variables are already within range, minimum is 176.8 and max 1732. Median 734.6
summary(data_long$creat_pre_mgdl) #variables are already within range, minimum is 2 max 19.572. Median 8.31
summary(data_long$creat_pre_mgl) #variables are already within range, minimum 20 max 195. Median 83

summary(data_long$urea_pre_gl) #range is from 0.1502 to 3.0 median 1.1652
summary(data_long$urea_pre_mmoll) #range is from minimum of 2.115 to max of 47.79 with median of 20
summary(data_long$urea_pre_mgdl) # range is from minimum 12.7 to max 287 median of 120. 

#it is possible that one unit (e.g. umoll) has a value and another unit has a missing value there
#Therefore, we fill creat micromol and we fill urea mgdL with any other unitvariables if they are not missing
data_long <- data_long %>%
    mutate(
        #fill creat
        creat_pre_umoll = case_when(
            !is.na(creat_pre_umoll) ~ creat_pre_umoll,
            is.na(creat_pre_umoll) & !is.na(creat_pre_mgl) ~ creat_pre_mgl  * 8.84017,
            is.na(creat_pre_umoll) & !is.na(creat_pre_mgdl) ~ creat_pre_mgdl * 88.4017,
        ),
        sex = case_when(
            sex == "Male" ~ 0, 
            sex == "Female" ~ 1,
            .default = NA),
        #sometimes a negative value is entered, like -99 which actually represents a missing value
        urea_pre_mmoll = if_else(urea_pre_mmoll < 0, NA, urea_pre_mmoll),
        ktv_pre = if_else(ktv_pre <= 0, NA, ktv_pre),
        weight_post = if_else(weight_post <= 0, NA, weight_post),
        #fill urea
        urea_pre_mgdl = case_when(
            !is.na(urea_pre_mgdl) ~ urea_pre_mgdl,
            is.na(urea_pre_mgdl) & !is.na(urea_pre_mmoll) ~ urea_pre_mmoll * 6.006,
            is.na(urea_pre_mgdl) & !is.na(urea_pre_gl) ~ urea_pre_gl*100),
    ) %>%
    group_by(id) %>%
    fill(age, sex, height, .direction = "downup") %>%
    ungroup()

summary(data_long$creat_pre_umoll) #NA 14507 so filling it with the other variables did not change 
summary(data_long$urea_pre_mgdl) #NA 12080 so 15 NA's were resolved by converting from other urea units. Median min and max have not changed.

#To calculate nPCR we need the BUN in stead of urea
# In the United States and a few other countries, the urea level in plasma or serum is measured as nitrogen and is called "Blood Urea Nitrogen" or BUN. The unit used for BUN is mg/dL.
# In the rest of the world, the whole urea molecule is measured, not just the nitrogen, and is reported in Standard International units of mmol/L.
# This measurement is about twice as high as the BUN measurement because BUN only measures the nitrogen part of the molecule (a Molecular Weight of 28), while urea measures the whole molecule (a Molecular weight of 60). Thus urea is approximately twice that of BUN (60/28 = 2.14)
#reference for conversion value of 2.14: "Clinical practice guidelines for nutrition in chronic renal failure. KDOQI, National Kidney Foundation. Am J Kidney Dis. 2000 Jun;35(6 Suppl 2):S17-S104. doi: 10.1053/ajkd.2000.v35.aajkd03517. Erratum in: Am J Kidney Dis 2001 Oct;38(4):917.
data_long <- data_long %>%
    mutate(
        BUN = urea_pre_mgdl / 2.1428)

#save dataset
saveRDS(data_long, file ="/Users/sroos6/Documents/Projecten/HDF pooling trials/data/convince_data_long_nutrition.rds")
