/>Explanation/guide for the coding for post-hoc analysis of CONVINCE
By Sanne Roos 
R files:

**-0-nutritionConvince-datacleaning-lab **

we perform the following adjustments:
1.	Make all variables lower case
2.	Make dataset from wide to long, where preposition “v01” or “scr” is split into variable visit and numbers from 0 to 17
3.	Creatinin has 3 seperate variables, one for each unit. Sometimes one of these units is missing but the other one is not. We try to complete creat micromol and if it is missing, we calculate it from the other units (mgl or mgdl) if these are not missing. We use following calculation factors: from mgL to micromol: multiply by 8.84017, from mgdL to micromol: multiply by 88.4017. 
4.	We do the same for urea: from mmolL to mgdL multiply by 6.006, from gLto mgdL multiply by 100. (Note: we convert to mgdL in stead of micromol which we used for creatinin. This makes calculating BUN easier).
5.	Some values are negative for example -99 in stead of NA, we fix this issue by replacing negative values by NA
6.	We make sure that age sex and height are filled for each row for each patient ID, as these do not change over the visits.
7.	We calculate a new variable called BUN (blood urea nitrogen) which we get from dividing urea in mgdL by 2.1428.
8.	(#reference for conversion value of 2.14: Tantray, Javeed & Mansoor, Sheikh & Choh Wani, Rasy & Nissa, Nighat. (2023). Estimation of blood urea. 10.1016/B978-0-443-19174-9.00009-X. 
9.	Finally, we save this to file called: convince_data_long_nutrition.rds

**-1-nutritionConvince-dataprocessing**
In which we use the following dataset:
 The one created in 0-nutritionConvince-datacleaning, which we saved under convince_data_long_nutrition.rds
And we perform the following adjustments:
1.	Calculate variable nPCR, SCI and LTI. Calculate baseline variables for each row (we need to recalculate these after imputation, but for inspection we create them in our dataset which includes the missing values)
2.	Datacleaning: 
•	Some dates are noted as “2999-01-01” or “2666-01-01”. Some of these visits still have data in that row. Replacing it as missing would be wrong. Therefore we use the following rules:
If both the previous and the next date are not missing, we replace it by the middle date in between those dates;
If only the previous date is available, we add 2 months (mean time in between the visits) from the previous visit date;
If both previous and next date are missing, we make the date missing (this happens one time, but this row has no data in it).
•	Some visit dates are before the randomisation date, futhermore, some visit dates are not in consecutieve order. We use the following rules: 
-First problem: visit 1 does not come after visit 0, both dates are after randomisation, and further follow-up dates are in order. 
#Solution1 : swap the rows of visit 1 and 0.
-Second: There are two dates which does not come after the other, it is not visit 1. (it seems a typo in the year, month is correct). 
#Solution2: add 2 months from previous date if this is not missing.
-Third: There are two cases in which the visit 0 and or 1 come before the randomisation date, these patients do not have further follow up as they withdrew from the study before it started.
#Solution 3: if dates are before randomisation in visit 0, make them equal to randomisation date.
-Fourth: There is one case where visit 0 en 1 are before randomisation date, but this case does have adequate folluw-up after that. 
#Solution 4: following solution 3, visit 0 now becomes the randomisation date. Visit 1 is then before visit 0, en then solution 1 comes and swaps visit 0 en 1.
3.	Make visit a time variable, calculate time from randomisation date (ran_date)
4.	Keep visits where visit date is not missing. (only keep rows for valid visits for that patient). 
5.	For nPCR, complete case analysis would result in dropping 435 patients, imputation is therefore better.
6.	Recalculate new Follow-up time and status variable (accounting for competing risks)
In the current database, status is death or alive at end of study. No competing risks. We want to censor if Ntx occurs or death or quitting dialysis, as measurments will not occur anymore afterwards. In our case, if people are transplanted, nPCR will not be measured anymore. They are "dropped" for our outcome. But transplantation also affects nPCR/nutrtional status/ultrafiltration rate (other publication). 
Status is 1 if death occurs, 0 if patient is censored/alive at study end without ntx/quitting dialysis, 2 if transplantation occurs/quittind dialysis as a competing risk (note that in the final analysis, we take 2 and 1 together, and do not seperately take competing risk as a seperate outcome).
Currently, follow-up is until last monitored, patients were monitored until april 2023 (end of study). Note that this is different compared to our analysis as in CONVINCE, people were still followed after NTx to follow up for mortality (but we censor them).
7.	We save this file to convince_nutrition_1.rds

**-2-Convince-imputation**
In which we use the following dataset:
 The one created in 1-nutritionConvince-dataprocessing, which we saved under convince_nutrition_1.rds"
Furthermore, for the prediction matrix we want the comorbidities to improve prediction for imputation. These are not filled in in this dataset (all NA). Therefore, we also load file CONVINCE_export_20230417.xlsx with comorbidity data.
And we perform the following:
-we select the rows of interest with comorbidity variables from CONVINCE_export_20230417, and left join these, matching by id.
-we filter on valid visit dates. We do not want imputation over visits that did no occur
Futhermore, only one patient has visit 16, we will restrict analyses to visit 15.
-we select variables of potential interest to predict missing variables. Note that we do not select calculated variables created in 1-nutritionConvince-dataprocessing (such as nPCR, UFR, LTI, SCI, BUN). It is better to recalculate them after imputating the variables that they consist of (such as kt/v etc).
-We fill dialysis vintage for each row, as this does not change over time and we don’t want that to get imputated. Furthermore, we set group to 1 (HDF) and 0 (HD), which is more easy to interpret for the imputation model).
-We added name to the imputation model, which is the name of the center. This contains spaces, therefore we remove spaces from the variable name. 
-We set hdf convective volume to 0 instead of NA for everyone on HD, otherwise the imputation model will estimate something for HD while it should all be 0. 
-We set up the predictor Matrix. We see that there are 95 ranks (number of linearly independent columns), and 202 columns, meaning there are columns that are linearly dependent on another. Looking at the logged Events in mice object, we see that name is constant. However, with and without name gives same warnings. Leaving some variables out do not change the warnings. After we inspect the graphic plot later on after the imputation, the results look good. We ignore the warning.
-We make vector nimp in which we define the variables that are predictive and that do not have to be imputated. We define dates. We define limp, the variables we want to impute (that go into the formula of nPCR, LTI, SCI, UFR). We set the variables that do not have to be imputated (including the dates) to 0 in the matrix. Cluster variable is set to -2, in our case we want the model to cluster for id.
-We set the method for imputation for the variables of our interest to 2l.pmm, and for the dates we specify we do not use any method, so mice “knows” dates do not have to be imputed (otherwise it gives error because mice cannot impute dates).
-We perform imputation with 10 datasets and 50 iterations.
After running the model with and without name in it, we get the same warnings of failed convergence. Thus, we keep name. 
NOTE: running the imputation results in several warnings regarding scaling, this does not effect the quality of imputation result, we ignore this.
We also get the message that model fails to converge in some of the iterations. Plotting the result shows random plots for imputation, means we are satisfied with the result.
-We save this imputation result to 
“imputation_object.Rdata”
-We turn it back into a dataframe and recalculate bmi, bun, npcr, sci, lti and ultrafiltration rate.
-We save this to “imputed_data_convince.Rdata”

**3-Convince-IPCW**
We load “imputed_data_convince.Rdata” which we created in 2-Convince-imputation
We only want to use the imputed datasets so we filter out the imputation 0 (with NA’s/missings in it).
We create a variable called inf_cens for informative censoring. The last row of each patient should be 1 if they stopped the study before the end for an informative reason. Stopping dialysis, kidney transplantation, clinical reasons, changed modality, patient decision or patient died are informative. (there are no more measurement points available and this is due to a reason that could have changed the nutritional value or UFR (other project) of interest). For example, patients with a worse nutritional status could die sooner and that could be the reason why they have no more npcr measurements. 
Just for the formality, we also create variable non_inf_cens which is for non informative censoring such as participant moved or Completed treatment.
We recalculate baseline values for our mixed model later on, and we recalculate time as a continuous variable (time was not in the imputation process so we recalculate it).
Next, we want to inspect how many people have informative censoring at each point. It seems that at visit 15, no one has informative censoring. At each visit point, the amount of informative censoring is below 100. 
We save what we have so far to "data_unweighted.Rdata"
-we try to make a model across all visits, that predicts change of not being censored due to an informative reason. We estimate this based on clinical predictors. We calculate inverse probability ipcw. We see that after weighting, density curves show improved overlap. See also supplemental figure 1 of the publication. 
As the frequency of informative censoring was low, we also fit a generalized linear model per imputation but over all visits, and we perform sensitiviy analyses with that. 

**4-CONVINCE-nutrtition unweight**
In which we use “data_unweighted.Rdata” created in 3-IPCW script.
Part 1: linear mixed models to see if nutritional variables differ over time between HD and HDF 
-We turn visit into a factor, to check if linear association is appropriate. Furthermore, each visit is about 3 months later for each patient +/- 1 or 2 weeks. Therefore, turn i tinto time in months, derived from the visit variable. 
-we make a list of the imputations to loop over
-we apply linear mixed models, using visit_time as a factor, with in the model the baseline nutritional variable, and we cluster for id. 
-we perform this on each imputation and pool from package mice. (we also pooled manually, this resulted in same outcome so we continue in the rest of the script with automatic pooling (function pool) from the mice package). 
-We plot the predicted values using ggpredict. 
Part 2: mortality association of the nutritional variables
#nutrition is not randomized, therefore we corrext in function cox_nutr for several confounding variables.
We apply both linear cox function (cox_nutr). But we also check for potential non-linear relationships with function cox_spline_nutr. we use normalised splines in stead of penalised splines (penalised splines somehow did not work with weights, normalised and penalised are similar). We let it choose knots automatically with 3 degrees of freedom. 
-note that we let the model center in the rubin_rule_cox_spline function, not in the cox function itself. As we want to center over the total pooled result and not per imputation (then each imputation could have a different center, we want one center for the plot). However, it is possible that that value does not exist exactly, which is why we use this bit of code:
#with approx function, you return a list of points which linearly interpolate given data points
    #x and y are the numeric vectorns giving the coordinaties of the points to be interpolated
    #xout is an optional set of numeric values specifying where interpolation is to take place (see also https://www.r-bloggers.com/2023/08/mastering-data-approximation-with-rs-approx-function/)
            center <- approx(x = data[[nutr_var]], y = data$b, xout = center_at)$y

 **4-CONVINCE-nutrtition weighted good**
 In which we use "data_weight_allvisits_imp.Rdata" created in 3-IPCW 
 -For the LMM, we want to use weights but lmer or lme4 only has precision weights not leveled weights: https://american-institutes-for-research.github.io/WeMix/articles/Weighted_Linear_Mixed_Effects_Models.html
 -we note that svylme and svyglm are also options.
 -we test all 3 to see if they give similar estimates. They do. However, svylme and svyglm do not seem to support time/visit as a categorical variable. So we choose WeMix. WeMix by default needs 2 weights, we only have one. So we set the other weight to 1 for everyone. Furthermore, the output is not tidy (no clear coefficient or standard deviation in the output). Therefore we have to do that ourselves with tidy output.
 -For the weighted cox model, we can still use the same model but we specify weights (cox can deal with a single weight) and robust standard errors
 

**5.2-CONVINCE-baseline-derivation.qmd**
First step is baseline derivation
-We use the unimputed dataset (so we filter .imp=0 from “imputed_data_convince.Rdata”). We filter visit 0 to retrieve baseline information.
-For baseline table, we selected some relevant aspects, and we will further refer to more complete baseline information in the original publication of the CONVINCE trial. For now, we selected
#general information:     "country", "sex", "age",
#vitals:    "sbp_pre", "dbp_pre", (blood pressure)
#dialysis information:     "dial_vintage" (how long they have been on dialysis), "ktv_pre" (Kt/V as measure of dialysis adequacy), "blood_flow" (in the vascular access), "uf_vol" (ultrafiltration volume),
#laboratory:    "creat_pre_umoll", "bun" (blood urea nitrogen), "crp_pre_mgdl", "hb_pre_mmoll", "urea_pre_mgdl", 
#variables of interest in analaysis nutrition:    "npcr", "lti", "sci", "bmi"
-we plot histograms of each variable to inspect the normality, and then specify the variables which are not normally distributed. These will be represented as median with IQR.
You can use render function to get the table into a word document





