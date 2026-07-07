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
8.	(#reference for conversion value of 2.14: "Clinical practice guidelines for nutrition in chronic renal failure. KDOQI, National Kidney Foundation. Am J Kidney Dis. 2000 Jun;35(6 Suppl 2):S17-S104. doi: 10.1053/ajkd.2000.v35.aajkd03517. Erratum in: Am J Kidney Dis 2001 Oct;38(4):917)
9.	Finally, we save this to file called: convince_data_long_nutrition.rds


