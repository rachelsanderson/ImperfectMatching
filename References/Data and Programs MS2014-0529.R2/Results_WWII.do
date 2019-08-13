clear
clear all
set more 1

cap log close

log using Results_WWII.log, replace

use MP_WII.dta, replace

***********************************
*generate macros for different sets of regressors
*********************************

global mom "divorced widow husbandaway marst_miss"
global kid "childageyears numkids  maxage minage length_name"
global county10 "sei_mean sei_sd p_urban p_widows children_singlemom poverty fem_lfp child_labor val_farm"
global countyd "CID2-CID73"
global match "datemiss"
global state "S2-S11"
global state_year "manwrat ageent labage contschl gen_total char_tot tot_edu_schools state_miss work_required limited_duration  monthlyallowfirstchild monthlyallowaddchild"
global cohort "yob1 yob2"
global cohortd "BY2-BY26"

***************************************************************************
*means for the entire sample
*produce a table with means and p-values of the differences with rejected
***************************************************************************

*how many observations and how many people?
codebook mpid
tab accepted

global det "year yob childageyears numkids  maxage minage length_name widow divorced husbandaway marst_miss datemiss"

*********************************************************************************************************
*summary stats

sum $det

mat T = J(12,3,.)
local j=1
foreach i of global det {
sum `i'
ttest `i', by(accepted)
mat T[`j',1] = r(mu_1)
mat T[`j',2] = r(mu_2)
mat T[`j',3] = r(p)
local j=`j'+1
}
matrix rownames T = $det
frmttable using AppendicS8cwwii.doc, statmat(T) varlabels replace ///
	ctitle("", Rejected, Accepted, p-value)

gen age1940=1940-yob

	sum childageyears age1940, d



***********************************************************************************
*Main results for treated v rejected
*use unique/best match for time being
*need to update using Bo's new programs
***********************************************************************************

global outcome "height weight bmi"

*use linear spec for continous variables
local doit="replace"
foreach y of global outcome { 
sum `y' if  ~accepted
local M=r(mean)
reg `y' accepted, cluster(fips)
outreg2 accepted using Table7, sdec(3) br se  `doit' excel ctitle(" `y' Unique no controls OLS") title("Linear with unique matches") addstat(Mean, `M')
reg `y' accepted $kid  $mom  $match $countyd $state_year $cohortd, cluster(fips)
outreg2 accepted using Table7, sdec(3) br se  append excel ctitle("`y'  Unique all controls OLS") addstat(Mean, `M')
local doit="append"
}

*unique  logit for dummy variables
global outcome2 "illiterate underweight obese blackwwii"

foreach y of global outcome2 { 
sum `y' if  ~accepted
local M=r(mean)
logit `y' accepted, cluster(fips)
outreg2 accepted using Table7, sdec(3) br se  append excel ctitle(" `y' Unique no controls Logit") addstat(Mean, `M')
logit `y' accepted $kid  $mom  $match $countyd $state_year $cohortd, cluster(fips)
outreg2 accepted using Table7, sdec(3) br se  append excel ctitle("`y'  Unique all controls Logit") addstat(Mean, `M')
local doit="append"
}

*censored regressions for education
gen cens=0
replace cens=-1 if illiterate
replace cens=1 if student
cnreg educ_yrs accepted, cluster(fips) censored(cens)
outreg2 accepted using Table7, sdec(3) br se  append excel ctitle("education censored all matches no controls")
cnreg educ_yrs accepted $kid  $mom  $match $countyd $state_year $cohort, cluster(fips) censored(cens)
outreg2 accepted using Table7, sdec(3) br se  append excel ctitle("education censored all matches all controls")

*graphs
*education graph
		proportion educ_yrs if accepted==1
		estimates store accepted
		proportion educ_yrs if accepted==0
		estimates store rejected
coefplot accepted rejected, vertical recast(bar) barwidth(0.2) fcolor(*.5) citop title("Figure 8a: Education distribution in WWII records") ytitle("Proportion") saving(Figure8a, replace)


*bmi graph
kdensity bmi if accepted==1, addplot(kdensity bmi if accepted==0) legend(label(1 "Accepted") label(2 "Rejected")) title("Figure 8b: Distribution of BMI in WWII records")  saving(Figure8b, replace) note("graph from unique matches only")

stop
