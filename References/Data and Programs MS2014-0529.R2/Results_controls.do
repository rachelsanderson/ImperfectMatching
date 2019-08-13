clear
clear all
clear mata
set more 1
set memory 500m
set matsize 5000
program drop _all
cap log close

log using Results_controls, replace
use MP_controls

*orphans
kdensity ageatdeath2 if controltype==5, addplot(kdensity ageatdeath2 if accepted==0 || kdensity ageatdeath2 if accepted) legend (label(1 "Orphans") label(2 "Rejected") label(3 "Accepted")) title("Figure 6a: Age at death. MP boys and Census orphans") note("Unique matches only. Orphans defined as children living in institutions") saving(Figure6a, replace)

*drop states where abandoned/divorced are eligible
drop if state=="Colorado" | state=="Minnesota" | state=="Ohio" | state=="Wisconsin" 

*other poor mothers
kdensity ageatdeath2 if group2=="Poor mothers" & ~widow, addplot(kdensity ageatdeath2 if accepted & ~widow || kdensity ageatdeath2 if ~accepted & ~widow) title("Figure 6b: Age at death. Children of poor non-widows") subtitle("In states where divorced/abandoned mothers are ineligible") legend(label(1 "Poor mothers") label(2 "Accepted") label(3 "Rejected")) note("Unique matches only. We exclude CO, MN, OH and WI")  saving(Figure6b, replace)

