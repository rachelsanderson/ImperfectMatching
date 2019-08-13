program MultiMatch_x,eclass
* version 11.2
/* 
	the convension for the variables will be :
		the first is a dependent variable
		the next k1 are x-variables measured in the y-data. They will be used in the reduced form
		the next k2 are the covariates in the equation of interest. They are typically measured in the x-data. They can overlap with the first k1
		the next k3 are the exogenous variables (incl instruments). This is only relevant for the 2sls
		the last is the identifier
	
	Options
	
		modeltype()
			0: OLS
			1: logit
			2: IV
			nothing else implemented
	
*/
if replay() {
		display "Replay not implemented"
	}
	else {
		syntax varlist(min=2) [if] [in] [, id(varlist), cluster(varlist), rf(varlist), iv(varlist), CLUSTERDummies, DETails, DEBug,  modeltype(real 0), ALLOWError,  Wls bs(real 0)]
	
		local nvars: word count  `varlist'
		local depvar: word 1 of `varlist'
		local regs: list varlist - depvar
		if !(`modeltype'==2) local iv=""
		local k1 : word count  `rf'
		local k2 : word count  `regs'
		local k3 : word count  `iv'
		local k0 : word count  `cluster'
		local k_id : word count `id'
		
		
		if !(`k_id'==1) {
			display " "
			display " "
			display " **** --- specify the identifier --- ****"
			display " "
			display " "
			exit
			}

		if !(`k0'==1) local cluster `id' 
		if ( (!(`modeltype'==1))&(`bs'==0)) local bs=200
		 
		
		if (`k1'>0){

			marksample touse

			local printdetails=1
			if !("`details'"=="details") local printdetails=0
			local debugging=1
			if !("`debug'"=="debug") local debugging=0
			local allowkappa=1
			if !("`allowerror'"=="allowerror") local allowkappa=0
			local gls=1
			if !("`wls'"=="wls") local gls=0
			local cludum=1
			if !("`clusterdummies'"=="clusterdummies") local cludum=0
			if ((`cludum'==1)&(	`k0'==0)) {
			display " "
			display " "
			display " **** --- clusterdummied requires clustering --- ****"
			display " "
			display " "
			exit
			}
			
			tempname b V nn ch2 kk ideleted
			

			mata: estimate_multi_match("`varlist'", "`rf'", "`iv'", "`id'", "`cluster'" ,"`touse'",`printdetails',`debugging',`k1',`k2',`k3', `modeltype', `allowkappa', `gls', `bs', `cludum')

			
			mat b = r(beta)
			mat V= r(V)
			matrix `nn'= r(N)
			matrix `kk'= r(K)
			matrix `ideleted' = r(I)
			local k1=`kk'[1,1]
			local k2=`kk'[1,2]
			local k3=`kk'[1,3]
			local k4=`kk'[1,4]
			local nx=`kk'[1,5]	
		
			display " "
			display as text _dup(78) "="
			display " "
			display "         Estimation Results using Procedure Described in ???:  "
			display " "
			if (`modeltype'==0) display "         Estimation of Linear Regression with Multiple Matches"
			if (`modeltype'==1) display "         Estimation of Logit with Multiple Matches"
			if (`modeltype'==2) display "         Estimation of IV Model with Multiple Matches"
			display "         Number of observations identified:  " %5.0f  `nn'[2,1] 
			display " "
		 
		 
			local  ddd="constRF"
			local i
			forvalues i=1/`k1' {
			local expvar: word `i' of `rf'
			local ddd  `ddd' `expvar'RF 
			}
			local ddd `ddd' "_cons" 
			local droppedv="" 
			forvalues i=1/`nx' {
			local ip1=`i'+1
			local expvar: word `ip1' of `varlist'
			if(`ideleted'[`i',1]>0.5){
			local droppedv  `droppedv' `expvar' 
			}
			else{
			local ddd  `ddd' `expvar'  
			}
			}
			
			if  ("`allowerror'"=="allowerror") local ddd="`ddd' kappa"
			matrix colnames b= `ddd' 
			matrix colnames V= `ddd'   
			matrix rownames V= `ddd' 
		
			local N = `nn'[2,1]

			ereturn clear	
			ereturn post b V, depname(`depvar') obs(`N') esample(`touse')
			ereturn local cmdline `"`0'"'
			ereturn local cmd "bop"
			display " "
			ereturn display 
			display " "
			if (`bs'>0) {
				display "Clustered at: `cluster'"
				}
				else {
				display "Analytic Standard Errors"
			}
				
			display " "
			display as text _dup(78) "="
			local ndr: word count  `droppedv'

			if (`ndr'>0.5){
				display " "
				display "Collinearity detected: "
				display " "
				local i
				forvalues i=1/`ndr' {
					local expvar: word `i' of `droppedv'
					display "  `expvar' dropped"
				}
				display " "
				display as text _dup(78) "="
			}

			display " "
			display " "
			}
		else {
			display " "
			display " "
			display " **** --- THE REDUCED FORM MUST HAVE EXPLANATORY VARIABLES IN THIS VERSION --- ****"
			display " "
			display " "
		}
	}
end
