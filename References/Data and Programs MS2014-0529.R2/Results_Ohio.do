
clear
clear all
clear mata
set more 1
set memory 500m
program drop _all
cap log close

log using Results_Ohio, replace

quietly do MultiMatch_x.do
quietly do matafile_MultiMatch_x.do

use MP_Ohio
global mom "divorced husbandaway marst_miss"
global kid "childageyears length_name sib2-sib8  maxage minage"
global match "datemiss"
global state_year "manwrat ageent labage contschl gen_total char_tot tot_edu_schools work_required limited_duration  monthlyallowfirstchild monthlyallowaddchild"
global state_year2 "manwrat ageent labage contschl gen_total char_tot tot_edu_schools  work_required limited_duration  monthlyallowfirstchild monthlyallowaddchild"
global countyO "totalexpendituresonrelief amountoutsidereliefcounty totalchildrenchildrenhomes"

tab fips, gen(OCD)
global countyd2 "OCD2-OCD17"

*Table S8
global det2 "year yob childageyears numkids  maxage minage length_name widow divorced husbandaway marst_miss datemiss"

preserve
keep if nmatcheseither==1
mat T = J(12,3,.)
local j=1
foreach i of global det2 {
sum `i'
ttest `i', by(accepted)
mat T[`j',1] = r(mu_1)
mat T[`j',2] = r(mu_2)
mat T[`j',3] = r(p)
local j=`j'+1
}
matrix rownames T = $det2
frmttable using tableS8.doc, statmat(T) varlabels replace ///
	ctitle("", Rejected, Accepted, p-value)	
restore


***********************************************
*graph differences 
kdensity ageatdeath2 if accepted==1 & nmatcheseither==1 & ageatdeath2>20, addplot(kdensity ageatdeath2 if accepted==0 & nmatcheseither==1 &  ageatdeath2>20)  legend(label(1 "Accepted") label(2 "Rejected") ) title("Figure 4: Distribution of age at death in Ohio") note("Unique matches only") saving(Figure4, replace)


******************************
*results Table S7a

*use only DMF matches

sum survivedto70a if ~accepted 
local M=r(mean)

MultiMatch_x survivedto70a accepted $kid  $mom $match  $countyd2 $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)  details 
outreg2 using TableS7b,  sdec(3) br se  excel keep(accepted) replace addstat(Mean, `M') ctitle("no county controls") 

MultiMatch_x survivedto70a accepted $kid  $mom $match  $countyO $countyd2 $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)   details
outreg2 using TableS7b,  sdec(3) br se  excel append keep(accepted)  ctitle("Add county controls")

sum survivedto70a if ~accepted & source_deathinfo=="Death Master File"

local M=r(mean)
MultiMatch_x survivedto70a accepted $kid  $mom $match $countyO $countyd2 $cohortd if source_deathinfo=="Death Master File". , modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)   details
outreg2 using TableS7b,  sdec(3) br se  excel append keep(accepted) addstat(Mean, `M') ctitle("Survived to `y' MP DOB county FE drop non-matches")

*repear with all matches

sum survivedto70 if ~accepted 
local M=r(mean)

MultiMatch_x survivedto70 accepted $kid  $mom $match  $countyd2 $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)  details
outreg2 using TableS7b,  sdec(3) br se  excel keep(accepted) `doit' addstat(Mean, `M') ctitle("no county controls") 

MultiMatch_x survivedto70 accepted $kid  $mom $match  $countyO $countyd2 $cohortd, modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4) cluster(fips)   details
outreg2 using TableS7b,  sdec(3) br se  excel append keep(accepted)  ctitle("Add county controls")

drop if ageatdeath2==. 
sum survivedto70 if ~accepted 
local M=r(mean)
MultiMatch_x survivedto70 accepted $kid  $mom $match $countyO $countyd2 $cohortd , modeltype(1)  id(mpid) rf(yob1 yob2 yob3 yob4)  cluster(fips)  details
outreg2 using TableS7b,  sdec(3) br se  excel append keep(accepted) addstat(Mean, `M') ctitle("Survived to `y' MP DOB county FE drop non-matches")



