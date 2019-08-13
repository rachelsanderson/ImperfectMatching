clear
cap log close

use "C:\Users\Anna_Aizer\Dropbox\Welfare\AER\forpublication\Data and Programs\iowafortable1b.dta",replace


/* table 1b*/

global t "iowatable1b.xls"
global cl "cluster(fips)"

reg found accepted, $cl 
outreg2 $o using $t, $o replace ti("Table 1b: Balance on Pre Treatment Characteristics" "Subsamples from Iowa and Ohio")
reg income accepted, $cl 
outreg2 using $t, $o append 
reg lninc accepted, $cl 
outreg2 $o using $t, $o append 
reg zeroinc accepted, $cl 
outreg2 $o using $t, $o append 
reg homeowner accepted, $cl 
outreg2 $o using $t, $o append 
reg homevalue accepted, $cl 
outreg2 $o using $t, $o append 
reg debt accepted, $cl 
outreg2 $o using $t, $o append 
reg yrsSchool accepted, $cl
outreg2 using $t, append $o 
reg literate accepted, $cl 
outreg2 $o using $t, $o append 
