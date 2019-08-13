/* this program produces table 1b Additional Retrospective Data on MP Applicants*/


clear
use "C:\Users\Anna_Aizer\Dropbox\Welfare\AER\forpublication\Data and Programs\ohiofortable1b.dta",replace

global t "census.xls"
global o "se bracket noaster excel label"

drop if censusyear==1930
global cl "cluster(fips)"

reg found accepted, $cl
outreg2 using $t,replace $o
reg native accepted, $cl
outreg2 using $t, append $o
reg owner accepted , $cl
outreg2 using $t, append $o
reg imputedincome accepted if owner~=., $cl
outreg2 using $t, append $o


keep if owner~=.
/* histograms*/
proportion imputedincome if accepted==1
        estimates store Accepted
        proportion imputedincome if accepted==0
        estimates store Rejected
        coefplot Accepted Rejected, vertical recast(bar) barwidth(1) fcolor(*.5) mlabel xlabel(none) citop xtitle("Imputed Income - Ohio") ytitle("Density") note(Accepted n=714 Rejected n=102) saving(ohio_h, replace)
