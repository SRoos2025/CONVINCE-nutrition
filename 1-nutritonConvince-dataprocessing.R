### datacleaning for dates and variables
### code by S. Roos
### last updated 12-18-2025

#0. set-up ----
#load packages
pacman::p_load ("rio",     #to recognize other files and load them in R
                "conflicted", #if a function is in more packages it makes you choose which one to use
                "tidyverse" #to work with nice syntaxes
)


#1.0 resolve package conflicts----
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

#1.1 read database, we use database created in script "0-nutritionConvince-datacleaning".----
data_long <- read_rds(file ="/Users/sroos6/Documents/Projecten/HDF pooling trials/data/convince_data_long_nutrition.rds")

#2.0 calculate nutritonal variables SCI LTI and npcr----
#calculate SCI (mg/kg/day) which we need for the LTI
#SCI (mg/kg/day) = 16.21 + 1.12 × [1 if male; 0 if female] – 0.06 × Age (years) – 0.08 × spKt/V urea + 0.009 × Pre-dialysis SCr (μmol/L)
#LTI Estimated LTI (kg/m²) = [(SCI × Post-HD Weight (kg) × 0.029) + 7.38] / [Height (m)2]
data_long <- data_long %>%
    mutate(
        #npcr works with urea BUN in mgdL, the formula below is for thrice weekly dialysis, we use midweek as we dont know exactly where the measurment took place in the week
        npcr = (BUN / (25.8 + (1.15 * ktv_pre) + 56.4/ktv_pre)) + 0.168,
        sci = (16.21 + if_else(sex == 0, 1.12, 0) - (0.06 *age) - (0.08 * as.numeric(ktv_pre)) + (0.009 * as.numeric(creat_pre_umoll))),
        lti = ((sci * weight_post*0.029)+7.38)/((height/100)^2),
        #make treatment factor for comparison later on
        group = factor(group, 
                           levels = c("High-flux hemodialysis", "High-dose hemodiafiltration"))
        )

#2.1 check missings----
#for each visit (before removing rows without valid visit date)
missings_before_row_selection <- data_long %>%
    group_by(visit) %>%
    summarise(percent_missing_npcr = mean(is.na(npcr)) * 100)

missings_before_row_selection <- data_long %>%
    group_by(visit) %>%
    summarise(percent_missing_sci = mean(is.na(sci)) * 100)

missings_before_row_selection <- data_long %>%
    group_by(visit) %>%
    summarise(percent_missing_lti = mean(is.na(lti)) * 100)

missings_before_row_selection <- data_long %>%
    group_by(visit) %>%
    summarise(percent_missing_lti = mean(is.na(bmi)) * 100)
# 
# missings_before_row_selection

#3.0 datacleaning----
#some rows have visit date "2999-01-01 or 2996-01-01". This happens 26 times: 
# sum(data_long$visit_date > as.Date("2025-01-01"), na.rm = TRUE)
#we want to inspect these dates to see if it is a typo or if these visit dates should all be NA (because all data is NA)
filter_wrong_visits <- data_long %>%
    filter(visit_date >= as.Date("2025-01-01"))
#inspecting this shows that for some patients there is still data present in rows with this date. 
#Simply replacing the visit date by NA would be wrong.

#so we make an inbetween date for those dates
data_long <- data_long %>%
    group_by(id) %>% #group by id because we will use lead and lag function below
    mutate(
        orig_date = as.Date(visit_date),# Keep original dates
        prev_date = lag(orig_date),
        next_date = lead(orig_date)
    )
#there is one person for wich previous and next date are both missing, but this row is empty:
# sum(is.na(filter_wrong_visits$prev_date) & is.na(filter_wrong_visits$next_date)) 

#apply rules to determine new date
data_long <- data_long %>%
    group_by(id) %>%
    mutate(
        visit_date = 
            case_when( 
                orig_date >= as.Date("2025-01-01") & #date is wrong
                    !is.na(prev_date) & !is.na(next_date) & #previous and next date not missing and both are correct 
                    prev_date <= as.Date ("2025-01-01") & next_date <= as.Date ("2025-01-01") ~  prev_date + (next_date - prev_date)/2, #if both previous and next visit date are available, take the middle.
                orig_date > as.Date("2025-01-01") & 
                    #if next date is missing or wrong, but previous date is not
                    !is.na(prev_date) & prev_date <= as.Date ("2025-01-01") & (is.na(next_date) | next_date >= as.Date("2025-01-01")) ~ prev_date + months(2),   # only previous visit date exists, add 2 months, this is the mean interval between visit dates
                #if both previous and next date are missing then NA because we can not determine at all when the visit took place.
                orig_date >= as.Date("2025-01-01") & 
                    (
                        (is.na(prev_date) & is.na(next_date)) |
                            (next_date >= as.Date("2025-01-01") | prev_date >= as.Date("2025-01-01"))
                    ) ~ NA,
              
                .default = as.Date(visit_date))
    ) %>%
    ungroup()

#There are now 0 dates >2025:
# filter_fixed_visits <- data_long %>%
# filter(visit_date > as.Date("2025-01-01")) %>%
# select(id, visit_date, orig_date, prev_date, next_date)

#you can inspect the new dates
# fixed_visits_inspect <- data_long %>%
#     filter(orig_date > as.Date("2025-01-01")) %>%
#     select(id, visit_date, visit, orig_date, prev_date, next_date)


#first create baseline variables for each row
data_long <- data_long %>%
    group_by(id) %>%
    mutate(
        #npcr
        baseline_npcr = npcr[visit == 0],
        #SCI
        baseline_sci = sci[visit == 0],
        #LTI
        baseline_lti = lti[visit == 0],
        #post dialysis weight
        baseline_weight_post = weight_post[visit == 0],
        #calculate time variable (in months) from visitdate
        visit_date = as.Date(visit_date),
        ran_date = as.Date(ran_date),
        end_eos_date = as.Date(end_eos_date),
        time = case_when(
            visit == 0 ~ 0,
            visit != 0 ~ as.numeric(ran_date %--% visit_date, unit = "months"))
    ) %>%
        group_by(id) %>%
            fill(baseline_npcr, baseline_sci, baseline_lti, baseline_weight_post, .direction = "downup") %>% #make sure baseline variables are present at each row
            ungroup() 

#this shows time can also be negative, so the visit is before randomisation date?
#summary(data_long$time)

#filter patients which have time from randomisation <0 (so visit date before randomisation date)
filter_wrong_time_after_randomisation <- data_long %>%
    filter(time < 0) %>%
    select(id, visit_date, visit, orig_date, prev_date, next_date, time, ran_date)

#get the patients that have at least one non-valid date, so visit date is not larger then previous date.
patients_wrong_visit_order <- data_long %>%
    arrange(id, visit) %>%
    group_by(id) %>%
    mutate(date_diff_ok = visit_date > lag(visit_date)) %>%
    summarise(valid = all(is.na(date_diff_ok) | date_diff_ok)) %>%
    filter(!valid) %>%  
    pull(id)
length(unique(patients_wrong_visit_order)) #7 cases



#look at the wrong ones, both the ones with visit date before randomisation date and the ones where visit order is not good.
#these are total 9 cases.
df_wrong <- data_long %>%
    arrange(id, visit) %>%
    group_by(id) %>%
    mutate(date_diff_ok = visit_date > lag(visit_date)) %>%
    ungroup() %>%
    filter((id %in% patients_wrong_visit_order | id %in% filter_wrong_time_after_randomisation$id)) %>%
    select(id, visit_date, visit, orig_date, prev_date, next_date, time, ran_date, date_diff_ok)

#What goes wrong in these  cases?
#1. visit 1 does not come after visit 0, but both dates are after randomisation, and further follow-up dates are in order. 
#Solution 1: swap visit 1 and 0


#2. There is are dates which does not come after the other, it is not visit 1. (it seems a typo in the year, month is correct). 

#Potential solution: add 2 months from previous date

#3. There are 3 cases in which the visit 0 and or 1 come before the randomisation date, some of these patients do not have further follow up as they withdrew. 
#Potential solution 3: if dates are before randomisation, make them equal to randomisation.

#4. There is one case where visit 0 en 1 are before randomisation date, but this case does have adequate folluw-up after that. 
#if you follow solution 3, visit 0 becomes randomisation date. Visit 1 is then equal to visit 0, en then solution 2 comes, adds two months. 

#5. There is one case where visit 2 == 3 and  (visit 7 = 8)
#Solution 2 fixes this.

#solution 3
data_long <- data_long %>%
    group_by(id) %>%
    mutate(
        visit0 = visit_date[visit == 0],
        visit1 = visit_date[visit == 1],
        visit_date = case_when(
            #solution 3, make dates equal to randomisation (we start with problem 3 so problem 4 solves itself with the solutions after)
            visit == 0 & visit_date < ran_date ~ ran_date,
            .default = visit_date)) %>%
    ungroup()

#now solution 1, we want to swap the entire rows not just the dates.
data_long <- data_long %>%
    group_by(id) %>%
    mutate(
        #create variable swap, to indicate when the rows need to be swapped
        swap = !is.na(visit0) & !is.na(visit1) & visit1 < visit0,
        #create another sorting key that can be integer 1 or 2
        sort_key = case_when(
            swap & visit == 0 ~ 2L,
            swap & visit == 1 ~ 1L,
            .default = visit)) %>%
            arrange(sort_key, .by_group = TRUE) %>%
    #now we can swap visit number
    mutate(
        visit = case_when(
           visit == 1 & !is.na(visit0) & !is.na(visit1) & visit1 < visit0 ~ 0,
           visit == 0 & !is.na(visit0) & !is.na(visit1) & visit1 < visit0 ~ 1,
           .default = visit))%>%
            select(-visit0, -visit1, -swap, -sort_key)%>%
            ungroup()

#now solution 2
data_long <- data_long %>%
    mutate(
        #recalculate previous visit date
        orig_date = as.Date(visit_date),               # Keep original dates
        prev_date = lag(orig_date),
        next_date = lead(orig_date),
        visit_date = case_when(
            #solution 2 add 2 months to previous date
            visit >= 1 & !is.na(prev_date) & visit_date <= prev_date ~ prev_date + months(2), #in this case we want 2 months and not the last daty of the nth month (roll back). So we don't use %m-%
            .default = visit_date)
    ) %>%
    ungroup()

#repeat rule 3 in case dates remain that are before randomisation date
data_long <- data_long %>%
    group_by(id) %>%
    mutate(
        visit_date = case_when(
            visit == 0 & visit_date < ran_date ~ ran_date,
            .default = visit_date)
    ) %>%
    ungroup()

# it is possible that ran date is now the same as the next visit date, so repeat rule 2
data_long <- data_long %>%
    mutate(
        #recalculate previous visit date
        orig_date = as.Date(visit_date),               # Keep original dates
        prev_date = lag(orig_date),
        next_date = lead(orig_date),
        visit_date = case_when(
            #solution 2 add 2 months to previous date
            visit >= 1 & !is.na(prev_date) & visit_date <= prev_date ~ prev_date + months(2), #in this case we want 2 months and not the last daty of the nth month (roll back). So we don't use %m-%
            .default = visit_date)
    ) %>%
    ungroup()

#calculate time again with correct dates
#there are 4 people with fup_y = 0 (this was already in loaded database)
#one person has one valid visit after visit 0, but it is less then a year before quitting
#so time should not be 0 there.
#for the other 3 with fup_y = 0, time should be 0 everywhere. After visit 0 all rows should be removed
data_long <- data_long %>%
    group_by(id) %>%
    mutate(
        time = case_when(
            visit == 0 ~ 0,
            fup_y == 0 & id != 186058 ~ 0, #to fix the 2 cases mentioned above
        visit != 0 & (fup_y != 0 | id == 186058) ~ as.numeric(ran_date %--% visit_date, unit = "months"))
    ) %>%
ungroup()

#remove the other visits of these 3 ID's but keep visit 0
data_long <- data_long %>%
    filter(!(id %in% c(203006, 263006, 190005) & visit != 0))

#inspect again
df_wrong_check <- data_long %>%
    arrange(id, visit) %>%
    group_by(id) %>%
    mutate(date_diff_ok = visit_date >= lag(visit_date)) %>% #this makes it easier to inspect if everything is correct, returns TRUE if everything went well.
    ungroup() %>%
    filter((id %in% patients_wrong_visit_order | id %in% filter_wrong_time_after_randomisation[["id"]])) %>%
    select(id, visit_date, visit, orig_date, prev_date, next_date, time, ran_date, date_diff_ok)

#everything is solved now
summary(data_long$time)



#4.0 keep rows with valid visit date, decide on imputation yes/no----
#keep rows where visit Date is not missing
data_long_valid_visits <- data_long %>%
    filter(!is.na(visit_date))

data_long_valid_visits <- data_long_valid_visits %>%
    mutate(
        npcr_missing = if_else(is.na(npcr), 1, 0)
    )

#925 have nPCR in all rows. So imputation is better
allrows <- data_long_valid_visits %>%
    group_by(id) %>%
    summarise(all_complete = all(npcr_missing == 0, na.rm = TRUE)) %>%
    filter(all_complete) %>%
    summarise(n_ids = n())

#5.0 recalculate FU time/status ----
#this shows that transplant flag works correctly, there are no missing flags when the reason of end of study is kidney transplantation
#test<- data_long %>%
    #filter(is.na(transplant_flag_eot) & end_eot_reason_recode_full == "Kidney transplant")

#make variable outcome, 2 for competing risk, 1 for death, 0 for censored
#note that afterwards, we combined 0 and 2 in the weighted cox model. 
data_long_valid_visits <- data_long_valid_visits %>%
    group_by(id) %>%
    mutate(
        transplant_flag_eot = if_else(!is.na(transplant_flag_eot) & transplant_flag_eot == "Yes", 1, 0), #convert flag to binary 1 or 0
        last_visit_date = max(visit_date),
        #if it is the last visit date and the reason treatment stopped was informative, set outcome to 2
        outcome = case_when(
            #end of treatment (EOT) reason
            visit_date == last_visit_date &
                end_eot_reason_recode_full %in% c(
                    "Stopped dialysis",
                    "Kidney transplant",
                    "Clinical reasons",
                    "Changed modality",
                    "Patient decision",
                    "Other"
                ) ~ 2,
            #if it is the last visit date and the reason treatment stopped was uninformative, set to 0
            visit_date == last_visit_date &
                end_eot_reason_recode_full %in% c(
                    "Completed treatment",
                    "Participant moved") &
                end_eos_reason_recode_full != "Patient died"~ 0, 
            # Death ONLY if there is not one of the EOT events below 
            #(completed treatment is OK, because this is filled in as the end_eot reason in case of death)
            visit_date == last_visit_date &
                !end_eot_reason_recode_full %in% c(
                    "Stopped dialysis",
                    "Kidney transplant",
                    "Clinical reasons",
                    "Changed modality",
                    "Patient decision",
                    "Other",
                    "Participant moved"
                ) &
                visit_date == last_visit_date & end_eos_reason_recode_full == "Patient died" ~ 1,
            .default = 0)
    ) %>%
    ungroup()

data_long_valid_visits <- data_long_valid_visits %>%
    group_by(id) %>%
    mutate(
        event_date = case_when(
            outcome == 1 ~ death_date,
            #if the outcome is censoring, there is no specific censordate so we take the last visit date
            outcome %in% c(0,2) ~ last_visit_date
        ),
        #event date is the last one, because in case of death outcome can be 0 or 2 before that
        event_date = max(event_date, na.rm = TRUE),
        time_to_outcome = as.numeric(ran_date %--% event_date, unit = "months")
    ) %>%
    ungroup()




#for cox model, we need intervals for observation, with outcome indicator at the end
data_long_valid_visits <- data_long_valid_visits %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(
        t_start = time,
        #t stop is the next visit time
        t_stop = lead(time),
        #make stop end at the last event (so where the next time is missing)
        t_stop = if_else(is.na(t_stop), time_to_outcome, t_stop),
        #event of only at last row number  
        endpt = case_when(
            #if it is the time to outcome and it is the last date, set the event to the corresponding number
            t_stop == time_to_outcome & visit_date == last_visit_date & outcome == 1~ 1, 
            t_stop == time_to_outcome & visit_date == last_visit_date & outcome == 2 ~ 2,
            t_stop == time_to_outcome & visit_date == last_visit_date & outcome == 0 ~ 0,
            .default = 0),
        status_end = 
            case_when(
                visit_date == last_visit_date ~ outcome),
        #fill this for each row
        status_end = max(status_end, na.rm=TRUE)
    ) %>%
    ungroup()

#if the endpt is censor/competing event (0/2), the eventdate is the last measurement date. 
#However, of course, it is not exactly at that date, it is somewhere after the last measurement but before the next visit.
#so on average, it is 1.5 month after the last measurement.
data_long_valid_visits <- data_long_valid_visits %>%
    group_by(id) %>%
    mutate(
        t_stop = case_when(
            visit_date == last_visit_date & t_start == t_stop ~ t_stop + 1.5,
            .default = t_stop)
    ) %>%
    ungroup()

#save file ----       
saveRDS(data_long_valid_visits, file ="/Users/sroos6/Documents/Projecten/HDF pooling trials/data/convince_nutrition_1.rds")
