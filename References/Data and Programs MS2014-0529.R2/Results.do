
clear
clear all
clear mata
set more 1
set memory 500m
set matsize 5000
program drop _all
cap log close

*********************************************************
*NOTES
*you need to install two programs in order to estimate the regressions that use
*all the matches to death certificates
******************************************
*NEW PROGRAM 

quietly do MultiMatch_x.do
quietly do matafile_MultiMatch_x.do

log using Results.log, replace

****************************************
**this file produces analysis comparing rejected and accepted applicants
**keeping only counties with both
**using stringent matches
**using new programs for OLS/IV with multiple matches
*****************************************

*note: there is more than one observation per person if the person was matched to more than one death certificate.
use MP_data, replace


***********************************
*generate macros for different sets of regressors
*********************************

global mom "divorced husbandaway marst_miss"
global kid "childageyears length_name sib2-sib8  maxage minage"
global county10 "sei_mean sei_sd p_urban p_widows children_singlemom poverty fem_lfp child_labor val_farm"
global countyd "CID2-CID75"
global match "datemiss"
global state "S2-S11"
global state_year "manwrat ageent labage contschl gen_total char_tot tot_edu_schools work_required limited_duration  monthlyallowfirstchild monthlyallowaddchild"
global state_year2 "manwrat ageent labage contschl gen_total char_tot tot_edu_schools  work_required limited_duration  monthlyallowfirstchild monthlyallowaddchild"
global cohort "yob1 yob2"
global cohortd "BY2-BY26"
global det "year yob childageyears datemiss numkids  maxage minage length_name widow divorced husbandaway marst_miss famearn miss nmatches keepssn ageatdeath2 logageatdeath"
global det2 "year yob childageyears datemiss numkids  maxage minage length_name widow divorced husbandaway marst_miss famearn miss nmatches"

***************************************************************************
*TABLE 2
*descriptive stats
*means and p-values of the differences with rejected
***************************************************************************

**FULL SAMPLE
*how many observations and how many people?
tab accepted
tab accepted if idtag

*how many families?
tab accepted if famidtag

*how many counties
codebook fips if idtag & accepted
codebook fips if idtag & ~accepted

*distribution by state (For Appendix Table SS)
ta state accepted if idtag
*how generous
sum amount_cpi if idtag==1, d
sum real_amount if idtag==1, d

*Sample with unique matches
tab accepted if nmatches==1
tab accepted if idtag & nmatches==1
tab accepted if famidtag & nmatches==1
codebook fips if idtag & nmatches==1 & accepted
codebook fips if idtag & nmatches==1 & ~accepted


*************************************
*compute mean differences at the individual level
*Table 1a reports the output 
************************************
preserve
*keep one observation per person
keep if idtag==1

*full sample
local doit="replace"

	foreach var in $det {
		reg `var' accepted, cluster(fips)  
		outreg2 using table1a_f.out, br se bdec(3)  `doit' 
local doit="append"
	
	}
	
*matched sample with age at death

local doit="replace"

	foreach var in $det {
		reg `var' accepted if nmatches==1, cluster(fips) 
		outreg2  using table1a_u.out, br se bdec(3)  `doit'
local doit="append"
	
	}


***********************************
*Appendix Figures S2
*invinvestigate differences in distributions 
************************************
			
		proportion childageyears if accepted==1
		estimates store accepted
		proportion childageyears if accepted==0
		estimates store rejected
		coefplot accepted rejected, vertical recast(bar) barwidth(0.2) fcolor(*.5) citop title("Figure S2a: Age Applied") ytitle("Proportion") saving(FigureS2a, replace)
		
		proportion numkids_b if accepted==1
		estimates store accepted
		proportion numkids_b if accepted==0
		estimates store rejected
		coefplot accepted rejected, vertical recast(bar) barwidth(0.2) fcolor(*.5) citop title("Figure S2b: Number of kids in family") ytitle("Proportion") saving(FigureS2b, replace)
					
		proportion minage_b if accepted==1
		estimates store accepted
		proportion minage_b if accepted==0
		estimates store rejected
		coefplot accepted rejected, vertical recast(bar) barwidth(0.2) fcolor(*.5) citop title("Figure S2c: Youngest kid in family") ytitle("Proportion") saving(FigureS2c, replace)
	

	gen log_fam=log(famearn)
label var log_fam "Log of predicted family income"
kdensity log_fam if accepted==1, addplot(kdensity log_fam if accepted==0) legend (label(1 "Accepted") label(2 "Rejected")) title("Figure S2d: Predicted income") note("Income predicted using the 1915 Iowa census. Incomes lower than or equal to 0 set equal to one.") saving(FigureS2d, replace)

	
*******************

restore
preserve

*Appendix table S7c
*check subsample representativeness. unique matches only
local doit replace
keep if nmatches==1
	foreach var in $det2 {
		reg `var' accepted if famearn<famearn_p50, cluster(fips) 
		outreg2 accepted using means_poor.out, br se bdec(3) nocons `doit' title("Below Median")
		reg `var' accepted if famearn>famearn_p50, cluster(fips) 
		outreg2 accepted using means_rich.out, br se bdec(3) nocons `doit' title("Above Median")

				reg `var' accepted if urban==1, cluster(fips) 
		outreg2 accepted using means_urban.out, br se bdec(3) nocons `doit' title("Urban")
		reg `var' accepted if urban==0, cluster(fips) 
		outreg2 accepted using means_rural.out, br se bdec(3) nocons `doit' title("Not Urban")

			reg `var' accepted if immig_county==1, cluster(fips) 
		outreg2 accepted using means_imm.out, br se bdec(3) nocons `doit' title("High p_imm")
		reg `var' accepted if immig_county==0, cluster(fips) 
		outreg2 accepted using means_native.out, br se bdec(3) nocons `doit' title("Low p_imm")

			reg `var' accepted if childageyears<=14, cluster(fips) 
		outreg2 accepted using means_less14.out, br se bdec(3) nocons `doit' title("less than 14")
		reg `var' accepted if childageyears<=10, cluster(fips) 
		outreg2 accepted using means_less10.out, br se bdec(3) nocons `doit' title("Less than 10")
		
		reg `var' accepted if numkids_b>=3 & numkids_b<=7, cluster(fips) 
		outreg2 accepted using means_size.out, br se bdec(3) nocons `doit' title("Family size 3-7")

			reg `var' accepted if yob>=1900 & yob<=1920, cluster(fips) 
		outreg2 accepted using means_old.out, br se bdec(3) nocons `doit' title("1900-1920")
		reg `var' accepted if yob>=1900 & yob<=1910, cluster(fips) 
		outreg2 accepted using means_oldest.out, br se bdec(3) nocons `doit' title("1900-1910")
				reg `var' accepted if yob>=1911 & yob<=1920, cluster(fips) 
		outreg2 accepted using means_youngest.out, br se bdec(3) nocons `doit' title("1911-1920")

		local doit append
	
	}

restore

	
********************************
*Table 3
*what do we know in this sample about reasons for rejection or end?
**********************************
ta reasonrejected2 if idtag
ta reasonends if idtag
	
*******************************************************************************
*does number of matches differ by accepted status?
******************************************************************************

tab manym accepted
tab nmatches if idtag==1, sum(accepted)
tab nmatches accepted if idtag==1

preserve
keep if idtag==1

*Figures 5a and 5b
*miss
sum miss if ~accepted
local M=r(mean)
logit  miss accepted  if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  replace excel ctitle("Basic no controls") title("Determinants of matching") addstat(Mean, `M')
estimates store Basic
logit  miss accepted $cohort $state_year if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("State and Cohort")
estimates store SC
logit  miss accepted $cohort $state_year $county10 if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("State, cohort and county 1910")
estimates store County
logit  miss accepted  $cohort $state_year $countyd if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("state, cohort and county dummies")
estimates store Countyd
logit  miss accepted $kid  $mom  $match $countyd $state_year $cohortd if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("Add individual controls")
estimates store Ind
coefplot Basic SC County Countyd Ind, yline(0) vertical keep(accepted) note("Graph displays coefficients and 95% confidence intervals.") title("Figure 5a: Does acceptance predict missing age at death?") subtitle("Coefficient from logit on Accepted=1") saving(Figure5a, replace) 

*manym
sum manym if ~accepted
local M=r(mean)
logit  manym accepted  if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  replace excel ctitle("Basic no controls") title("Determinants of matching") addstat(Mean, `M')
estimates store Basic
logit manym accepted $cohort $state_year if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("State and Cohort")
estimates store SC
logit manym accepted $cohort $state_year $county10 if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("State, cohort and county 1910")
estimates store County
logit manym accepted  $cohort $state_year $countyd if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("state, cohort and county dummies")
estimates store Countyd
logit manym accepted $kid  $mom  $match $countyd $state_year $cohortd if idtag==1, cluster(fips)
outreg2 accepted using selection,  sdec(3) br se  append excel ctitle("Add individual controls")
estimates store Ind
coefplot Basic SC County Countyd Ind, yline(0) vertical keep(accepted) note("Graph displays coefficients and 95% confidence intervals.") title("Figure 5b: Does acceptance predict multiple matches") subtitle("Coefficient from logit on Accepted=1") saving(Figure5b, replace) 

restore

****************************************************************
*Table 2
*determinants of who is accepted, who gets more money and how long they get money for.
***************************************************************

*Predicting acceptance and other generosity measures
global out2 "accepted logamount logdur"
sum $out2
local doit="replace"
foreach y of global out2 {
sum `y' if idtag==1 
local M=r(mean)
areg `y' $kid  $mom  $match $state_year $cohortd if idtag==1, absorb(fips) cluster(fips)
outreg2 childageyears sib2-sib8  maxage minage length_name divorced husbandaway marst_miss datemiss using table2,  excel sdec(3) br se  `doit' addstat(Mean, `M') ctitle(full `y') title("Determinants of generosity--full controls")
sum `y' if nmatches==1
local M=r(mean)
areg `y' $kid  $mom  $match $state_year $cohortd if nmatches==1, absorb(fips) cluster(fips)
outreg2 childageyears sib2-sib8  maxage minage length_name divorced husbandaway marst_miss datemiss using table2,  excel sdec(3) br se  append addstat(Mean, `M') ctitle(matched `y')

local doit="append"
}

***********************************************************************************
***********************************************************************************
*Main results for treated v rejected
***********************************************************************************

***************
***graphs******

*Figure 1
kdensity ageatdeath2 if accepted==1 & nmatches==1 & ageatdeath2>20, addplot(kdensity ageatdeath2 if accepted==0 & nmatches==1 & ageatdeath2>20) legend (label(1 "Accepted non-widows") label(2 "Rejected widows")) title("Figure 1: Distribution of age at death") subtitle("unique matches only") note("deaths below 20 dropped") saving(Figure1, replace)

*stochastic dominance test
ranksum ageatdeath2 if nmatches==1, by(accepted) porder

*Figure 3
*by income
kdensity ageatdeath2 if accepted==1 & nmatches==1 & ageatdeath2>20 & famearn<=famearn_p50, addplot(kdensity ageatdeath2 if accepted==0 & nmatches==1 & ageatdeath2>20  & famearn<=famearn_p50) legend (label(1 "Accepted") label(2 "Rejected")) title("Figure 3a: Predicted income below median in sample") note("unique matches only")  saving(Figure3a1, replace)

kdensity ageatdeath2 if accepted==1 & nmatches==1 & ageatdeath2>20 & famearn>famearn_p50, addplot(kdensity ageatdeath2 if accepted==0 & nmatches==1 & ageatdeath2>20  & famearn>famearn_p50) legend (label(1 "Accepted") label(2 "Rejected")) title("Figure 3b: Predicted income above median in sample") note("unique matches only")  saving(Figure3a2, replace)
*by rural status
kdensity ageatdeath2 if accepted==1 & nmatches==1 & ageatdeath2>20 & urban==1, addplot(kdensity ageatdeath2 if accepted==0 & nmatches==1 & ageatdeath2>20  & urban==1) legend (label(1 "Accepted") label(2 "Rejected")) title("Figure 3c: % urban above median in sample") note("unique matches only")  saving(Figure3b1, replace)
kdensity ageatdeath2 if accepted==1 & nmatches==1 & ageatdeath2>20 & urban==0, addplot(kdensity ageatdeath2 if accepted==0 & nmatches==1 & ageatdeath2>20  & urban==0) legend (label(1 "Accepted") label(2 "Rejected")) title("Figure 3d: % urban below median in sample") note("unique matches only")  saving(Figure3b2, replace)

*income

**********************
*MAIN RESULTS TABLE 4
*********************

*panel A: unique matches
preserve
keep if nmatches==1 

sum ageatdeath2 if ~accepted 
local M=r(mean)

reg logageatdeath accepted $state $cohortd, cluster(fips) 
outreg2  using Table4,  br se bdec(4) replace excel keep(accepted) addstat(Mean, `M') ctitle("Age at death state and cohort FE") title("Age at death")

reg logageatdeath accepted $kid  $mom  $match $county10 $state_year $cohortd, cluster(fips) 
outreg2  using Table4,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("Age at death all controls no county FE")

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd, cluster(fips) 
outreg2  using Table4,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("Age at death county FE")

reg logageatdeath_sa accepted $kid  $mom  $match $countyd $state_year $cohortd, cluster(fips) 
outreg2  using Table4, br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("Age at death SSA DOB")


*panel b: survival use all matches
restore
sum survived* accepted $kid  $mom $match $state_year nmatches

global ages2 "60 70 80"
*estimates for all

*estimates for all
foreach y of global ages2 {
sum survivedto`y'a if ~accepted
local M=r(mean)
MultiMatch_x survivedto`y'a accepted $state $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)  details
outreg2 using Table4,  sdec(3) br se excel append addstat(Mean, `M') ctitle("Survived to `y' DMF only state and cohort") keep(accepted)

MultiMatch_x survivedto`y'a accepted $kid  $mom $match $county10 $state_year $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)  details 
outreg2 using Table4,  sdec(3) br se excel append addstat(Mean, `M') ctitle("Survived to `y' full no C FE" ) keep(accepted)

MultiMatch_x survivedto`y'a accepted $kid  $mom $match $state_year $countyd $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)  details
outreg2 using Table4,  sdec(3) br se  excel append addstat(Mean, `M') ctitle("Survived to `y' full county and cohort FE plus all controls" ) keep(accepted)

MultiMatch_x survivedto`y' accepted $kid  $mom $match $state_year $countyd $cohortd , modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)  details
outreg2 using Table4,  sdec(3) br se excel  append addstat(Mean, `M') ctitle("Survived to `y' SSN DOB full county and cohort FE plus all controls") keep(accepted)

*unique matches only                       
sum survivedto`y'a if ~accepted & nmatches==1
local M=r(mean)
logit survivedto`y'a accepted $kid  $mom $match $state_year $cohortd $countyd if nmatches==1, cluster(fips)
outreg2 using Table4,  sdec(3) br se excel append addstat(Mean, `M') ctitle("Survived to `y' unique matches only") keep(accepted)


local doit="append"
}


*************************************************************************************************
*Table 5
***************************

**by income
keep if nmatches==1 
sum ageatdeath2 if ~accepted & famearn<famearn_p50
local M=r(mean)
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if famearn<famearn_p50, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) replace excel keep(accepted) addstat(Mean, `M') ctitle("log age at death poor")

sum ageatdeath2 if ~accepted & famearn>=famearn_p50
local M=r(mean)
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if famearn>=famearn_p50, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("log age at death rich")

*test equality
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if famearn<famearn_p50
estimates store Poor
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if famearn>=famearn_p50
estimates store Rich 
suest Poor Rich, cluster(fips)
test [Poor_mean]accepted=[Rich_mean]accepted

**by urban
sum ageatdeath2 if ~accepted & urban==1
local M=r(mean)
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if urban==1, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("log age at death urban")

sum ageatdeath2 if ~accepted & urban==0
local M=r(mean)
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if urban==0, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("log age at death rural")

*test equality
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if urban==1 
estimates store urban
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if urban==0
estimates store rural
suest urban rural, cluster(fips)
test [urban_mean]accepted=[rural_mean]accepted


**by share immigrants
sum ageatdeath2 if ~accepted & immig_county==1
local M=r(mean)
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if immig_county==1, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("log age at death immigra c")

sum ageatdeath2 if ~accepted & immig_county==0
local M=r(mean)
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if immig_county==0, cluster(fips) 
outreg2 accepted using Table5 ,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("log age at death not immigrant")

*test equality
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if immig_county==1
estimates store immi
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if immig_county==0
estimates store nativ

suest immi nativ, cluster(fips)
test [immi_mean]accepted=[nativ_mean]accepted

*test that other subgroups are the samea s main analysis

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd 
estimates store Main

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if childageyears<=14, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("<14")

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if childageyears<=10, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("<10")

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if yob<=1920, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("<=1920")

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if yob>=1900 & yob<=1910, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("yob>=1900 & yob<=1910")

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if yob>=1911 & yob<=1920, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("yob>=1911 & yob<=1920")

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if numkids>=3 & numkids<=7, cluster(fips) 
outreg2 accepted using Table5,  br se bdec(4) append excel keep(accepted) addstat(Mean, `M') ctitle("numkids>=3 & numkids<=7")

*test coefficients
reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if childageyears<=14
estimates store age14
suest age14 Main
test [Main_mean]accepted=[age14_mean]accepted

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if childageyears<=10
estimates store age10
suest age10 Main
test [Main_mean]accepted=[age10_mean]accepted

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if yob>=1900 & yob<=1920
estimates store old
suest old Main
test [Main_mean]accepted=[old_mean]accepted

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if yob>=1900 & yob<=1910
estimates store oldest
suest oldest Main
test [Main_mean]accepted=[oldest_mean]accepted

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if yob>=1911 & yob<=1920
estimates store youngest
suest youngest Main
test [Main_mean]accepted=[youngest_mean]accepted

reg logageatdeath accepted $kid  $mom  $match $countyd $state_year $cohortd if numkids>=3 & numkids<=7
estimates store fam
suest fam Main
test [Main_mean]accepted=[fam_mean]accepted

******************************
*Figure 2
*******************************

gen full_beta=.
gen full_se=.
gen full_rejectedmean=.

gen unique_beta=.
gen unique_se=.
gen unique_rejectedmean=.


local doit="replace"
forvalues y=58(1)84 {
gen sto`y'=(ageatdeath2>=`y')
replace sto`y'=0 if ageatdeath2==.

*full sample

sum sto`y' if ~accepted
local M=r(mean)
replace full_rejectedmean=`M'
MultiMatch_x sto`y' accepted $kid  $mom $match $state_year $countyd $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4)   details cluster(fips)
outreg2 using byage,  br se bdec(4) `doit' keep(accepted) excel ctitle("`y'") addstat(Mean, `M') title("Missing assumed dead") 
replace full_beta=_b[accepted] if ageatdeath2==`y'
replace full_se=_se[accepted] if ageatdeath2==`y'


*unique matches
sum sto`y' if ~accepted & nmatches==1
local M=r(mean)
replace unique_rejectedmean=`M'
logit sto`y' accepted $kid  $mom $match $state_year $countyd $cohortd if nmatches==1 , cluster(fips)
outreg2 using byageunique,  br se bdec(4) `doit' keep(accepted) excel ctitle("`y'") addstat(Mean, `M') title("Logit Unique matches") 
replace unique_beta=_b[accepted] if ageatdeath2==`y'
replace unique_se=_se[accepted] if ageatdeath2==`y'

local doit="append"
}
*save data for graph
preserve
keep unique_beta unique_se full_beta full_se full_rejectedmean unique_rejectedmean ageatdeath2
rename ageatdeath2 age
collapse unique_beta unique_se full_rejectedmean unique_rejectedmean full_beta full_se, by(age)
sort age
save coefficients_byage.dta, replace

restore

*repeat to do joint test
*se clustered in SUEST

*full sample
local doit="replace"
forvalues y=58(1)84 {
logit sto`y' accepted $kid  $mom $match $state_year $countyd $cohortd if nmatches==1 
estimates store age`y'

}
suest age58 age59 age60 age61 age62 age63 age64 age65 age66 age67 age68 age69 age70 age71 age72 age73 age74 age75 age76 age77 age78 age79 age80 age81 age82 age83 age84, cluster(fips)
test accepted


*****************************************************************************************
*TABLE S7a: Robustness checks
*do with 70 only
*****************************************************************************************
gen d=1
bysort mpid: gen counter=sum(d)

*unique matches
sum survivedto70a if nmatches==1 & ~accepted
local M=r(mean)
logit  survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd if nmatches==1, cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  replace excel ctitle("Unique matches only") addstat(Mean, `M')
*unique matches assuming nonmatches are dead (also counted as unique)
sum survivedto70a if (nmatches==1 | nmatches==0) & ~accepted
local M=r(mean)
logit  survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd if nmatches==1 | nmatches==0, cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  append excel ctitle("Unique matches only, missing counted as dead") addstat(Mean, `M')
*random match
sum survivedto70a if counter==1 & ~accepted
local M=r(mean)
logit  survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd if counter==1, cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  append excel ctitle("Random match") addstat(Mean, `M')
*all matches
sum survivedto70a if ~accepted
local M=r(mean)
logit  survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd, cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  append excel ctitle("All matches treated as observations") addstat(Mean, `M')
*highet quality match

sum survivedto70a if ~accepted & keepssn==1
local M=r(mean)
logit  survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd if keepssn==1, cluster(fips) 

*our method
sum survivedto70a if ~accepted 
local M=r(mean)
MultiMatch_x survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd, modeltype(1) id(mpid) rf(yob1 yob2 yob3 yob4) details cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  append excel ctitle("MLE using all matches") addstat(Mean, `M')

*matched on propensity score
preserve
codebook famid
gen childage=floor(childageyears)
replace childage=0 if childage<0
replace yob=floor(yob)
tab childage if counter==1, g(A)

local a=0
forvalues x=1(1)19 {
bysort famid: egen n_kids_age`a'=sum(A`x')
label var n_kids_age`a' "# kids in family of age `a'"
local a=`a'+1
}

xi: pscore accepted i.yob A* n_kids_age*  $mom $countyd $match $state_year $cohortd if counter==1, pscore(pscore) blockid(block) comsup
bysort mpid: egen pscore_i=mean(pscore)
gen matched=(pscore~=. | pscore==0)
drop pscore
rename pscore_i pscore
gen a_ps=accepted*pscore
bysort matched: sum $kid  $mom $county10 $match $state $cohort

	*keep matched sample and re-estimate
keep if matched
codebook fips
sum survivedto70a if ~accepted 
local M=r(mean)
MultiMatch_x survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd, modeltype(1) id(mpid) rf(yob1 yob2 yob3 yob4) details cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  append excel title("MLE with PS matched sample") addstat(Mean, `M')

restore

********************
*by number of matches
preserve
drop if nmatches>=3
sum survivedto70a if ~accepted
local M=r(mean)
MultiMatch_x survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd, modeltype(1) id(mpid) rf(yob1 yob2 yob3 yob4) details cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  append excel title("MLE dropping 3+ matches") addstat(Mean, `M')

drop if nmatches==0
sum survivedto70a if ~accepted
local M=r(mean)
MultiMatch_x survivedto70a accepted $kid  $mom  $match $countyd $state_year $cohortd, modeltype(1) id(mpid) rf(yob1 yob2 yob3 yob4) details cluster(fips)
outreg2 accepted using TableS7,  sdec(3) br se  append excel title("MLE only 1 or 2 matches") addstat(Mean, `M')

restore

*end robustness checks with full sample
*******************************************
preserve

**********************************************************************************************************************
*appendix FigureS1--missing data patterns
keep if idtag
collapse miss, by(yob)
label var yob "Year of birth"
label var miss "% missing age at death"
scatter miss yob, c(l) saving(FigureS1b, replace) title("Figure S1b: Pattern of missing data by cohort") subtitle("Among MP boys. Full sample")
restore

*******************************************************************************************************************************
*Table S7b
*AGGREGATE AT COUNTY LEVEL*YEAR AND RE-RUN

*with and without controls like in Table 6
*with and without pop weights
keep if (nmatches==1 | nmatches==0)

*collapse (mean) survivedto70a logageatdeath ageatdeath2 accepted $kid  $mom $match $state_year $state $countyd $cohortd (sum) nobs=d, by(fips year yob)
collapse (mean) survivedto70a logageatdeath ageatdeath2 accepted $kid  $mom $match $state_year $state $countyd $cohortd yob (sum) nobs=d, by(fips year)
codebook survivedto70a logageatdeath accepted yob nobs fips

*state-of-birth and YOB dummies only
sum survivedto70a if ~accepted
local M=r(mean)
reg survivedto70a  accepted $state $cohortd, cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  replace excel ctitle("Survived to 70 state-of-birth and YOB dummies") addstat(Mean, `M')

sum ageatdeath2 if ~accepted
local M=r(mean)
reg ageatdeath2  accepted $state $cohortd, cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  append excel ctitle("Age at death state-of-birth and YOB dummies") addstat(Mean, `M')

*all controls
sum survivedto70a if ~accepted
local M=r(mean)
reg survivedto70a  accepted $kid  $mom $match $state_year $countyd $cohortd, cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  append excel ctitle("Survived to 70 all controls") addstat(Mean, `M')

sum ageatdeath2 if ~accepted
local M=r(mean)
reg ageatdeath2  accepted $kid  $mom $match $state_year $countyd $cohortd, cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  append excel ctitle("Age at death all controls") addstat(Mean, `M')

*use number of observations as weights
*state-of-birth and YOB dummies only
sum survivedto70a if ~accepted [w=nobs]
local M=r(mean)
reg survivedto70a  accepted $state $cohortd [w=nobs], cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  append excel ctitle("Survived to 70 state-of-birth and YOB dummies WEIGHTS") addstat(Mean, `M')

sum ageatdeath2 if ~accepted [w=nobs]
local M=r(mean)
reg ageatdeath2  accepted $state $cohortd [w=nobs], cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  append excel ctitle("Age at death state-of-birth and YOB dummies") addstat(Mean, `M')

*all controls
sum survivedto70a if ~accepted [w=nobs]
local M=r(mean)
reg survivedto70a  accepted $kid  $mom $match $state_year $countyd $cohortd [w=nobs], cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  append excel ctitle("Survived to 70 all controls") addstat(Mean, `M')

sum ageatdeath2 if ~accepted [w=nobs]
local M=r(mean)
reg ageatdeath2  accepted $kid  $mom $match $state_year $countyd $cohortd [w=nobs], cluster(fips)
outreg2 accepted using TableS7b,  sdec(3) br se  append excel ctitle("Age at death all controls") addstat(Mean, `M')

*Figure2
clear
use coefficients_byage.dta, replace
sort age
drop if age<=57 
drop if age>=85

generate full_low = full_beta - 1.96* full_se
generate full_hi = full_beta + 1.96* full_se
gen full_me= full_beta*(1- full_rejectedmean)

label var full_beta "Beta from Logit"
label var full_me "Marginal effect as a fraction of rejected mean"
label var age "Age at Death"
label var full_low "low"
label var full_hi "hi"

twoway rarea full_low full_hi age , sort color(gs14)  || scatter  full_beta full_me age, c(l l) title("Figure 2a: Probability of surviving past given age") subtitle("Full sample, missing age at death assumed dead")
graph save Figure2a, replace

generate unique_low = unique_beta - 1.96* unique_se
generate unique_hi = unique_beta + 1.96* unique_se
gen unique_me= unique_beta*(1- unique_rejectedmean)

label var unique_low "low"
label var unique_hi "hi"
label var unique_beta "Beta from Logit"
label var unique_me "Marginal effect as a fraction of rejected mean"

twoway rarea unique_low unique_hi age , sort color(gs14) || scatter  unique_beta unique_me age, c(l l) title("Figure 2b: Probability of surviving past given age") subtitle("Matched sample. Unique matches only")
graph save Figure2b, replace


