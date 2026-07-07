/>Explanation/guide for the coding for post-hoc analysis of CONVINCE
By Sanne Roos 
R files:
-0-nutritionConvince-datacleaning-lab 

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

-1-nutritionConvince-dataprocessing
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



