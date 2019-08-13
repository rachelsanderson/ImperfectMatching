clear
clear all
clear mata
set more 1
set memory 500m
program drop _all
cap log close

log using results_1940.log, replace
********************************************************************************************
**this file produces analysis comparing rejected and accepted applicants
**using stringent matches. only one match per person with preliminary release of 1940 census.
**********************************************************************************************

use MP_1940, clear
sum

***********************************
*generate macros for different sets of regressors
*********************************

global mom "divorced husbandaway marst_miss"
global kid "childageyears length_name sib2-sib8  maxage minage"
global county10 "sei_mean sei_sd p_urban p_widows children_singlemom poverty fem_lfp child_labor val_farm"
global countyd "CID2-CID46"
global match "datemiss"
global state "S2-S11"
global state_year "manwrat ageent labage contschl gen_total char_tot tot_edu_schools state_miss  work_required limited_duration  monthlyallowfirstchild monthlyallowaddchild"
global cohort "yob1 yob2"
global cohortd "BY2-BY26"

sum $kid  $mom  $match $countyd $state_year $cohortd income educ black

***************************************************************************
*means for the entire sample
*produce a table with means and p-values of the differences with rejected
***************************************************************************

*how many observations and how many people?
tab accepted

*compute mean differences at the individual level. Table S8
global det "year yob childageyears numkids  maxage minage length_name widow divorced husbandaway marst_miss datemiss"

sum $det

mat T = J(12,2,.)
local j=1
foreach i of global det {
sum `i'
qui ttest `i', by(accepted)
mat T[`j',1] = r(mu_1)
mat T[`j',2] = r(mu_2)
local j=`j'+1
}
matrix rownames T = $det
frmttable using AppendixS8c.doc, statmat(T) varlabels replace ///
	ctitle("", Rejected, Accepted)
	

***********************************************************************************
*Main results for treated v rejected
***********************************************************************************

gen log_inc=log(income)
label var log_inc "Log of Annual Income in 1939"
***************
***graphs******
kdensity log_inc if accepted==1, addplot(kdensity log_inc if accepted==0) legend (label(1 "Accepted") label(2 "Rejected")) title("Distribution of log income 1940 census")  saving(Figure7a, replace)

		proportion educ if accepted==1
		estimates store accepted
		proportion educ if accepted==0
		estimates store rejected
coefplot accepted rejected, vertical recast(bar) barwidth(0.2) fcolor(*.5) citop xtitle("Education distribution in 1940 census") ytitle("Proportion") saving(Figure7b, replace)

**********************
*MAIN RESULTS TABLE 6
*********************
global outcome "income educ black"
*use  linear speac

local doit="replace"
foreach y of global outcome { 
sum `y' if  ~accepted
local M=r(mean)
reg `y' accepted, cluster(fips)
outreg2 accepted using Table6, sdec(3) br se  `doit' excel ctitle(" `y' no controls OLS") title("Linear with unique matches") addstat(Mean, `M')
reg `y' accepted $kid  $mom  $match $countyd $state_year $cohortd, cluster(fips)
outreg2 accepted using Table6, sdec(3) br se  append excel ctitle("`y' all controls OLS") addstat(Mean, `M')
local doit="append"
}

sum income, d
sum income if ~accepted, d
local p25=r(p25)
local p50=r(p50)
local p75=r(p75)

gen above25=(income>`p25')
gen above50=(income>`p50')
gen above75=(income>`p75')
replace above25=. if income==.
replace above50=. if income==.
replace above75=. if income==.

sum above*

		logit above25 accepted ,  cluster(fips) 
		estimates store p25
		logit above50 accepted ,  cluster(fips) 
		estimates store p50
		logit above75 accepted ,  cluster(fips) 
		estimates store p75
		logit above25 accepted $kid  $mom  $match $countyd $state_year $cohortd,  cluster(fips) 
		estimates store pc25
		logit above50 accepted $kid  $mom  $match $countyd $state_year $cohortd,  cluster(fips) 
		estimates store pc50
		logit above75 accepted $kid  $mom  $match $countyd $state_year $cohortd,  cluster(fips) 
		estimates store pc75
	
coefplot p25 pc25 p50 pc50  p75 pc75, yline(0) vertical keep(accepted) note("Graph displays coefficients and 95% confidence intervals, with and without controls") title("Fifgure S5: Coefficient of Accepted") subtitle("on probability of having income above percentile") saving(FigureS5, replace)  

