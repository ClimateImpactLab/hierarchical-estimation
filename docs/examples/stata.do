/*
Estimates global mortality response by year*adm2*age

:author: My Name
:contact: `myname@gmail.com <mailto:myname@gmail.com>`_
:modified: 2016/08/19
:team: Mortality
:lead: Amir Jina


Input variables
~~~~~~~~~~~~~~~

global_mortality_admin2_year_3agegrp.dta
    :location: ~/Dropbox/ChicagoRAs_Reanalysis/interpolation/data/consolidated
    :description: Global baseline mortality data
    :version: GCP.Mortality.GlobalBaselineMortality.2016-08-19


Output variables
~~~~~~~~~~~~~~~~

global_6countries_response_adm0prcp_global_year_FE.dta
    :location: ~/Dropbox/GCP/MORTALITY/tables/Table_2/Panel_A/ster
    :description: 6-country response surface with year FE
    :version: ``GCP.Mortality.Global_6countries_response_adm0prcp_yearFE.2016-08-19``

    The year FE dataset uses a linear regression with global year fixed effects:

        :math:`M_{i, t}=\sum_{i\in K}{\hat{\beta}_{i, t}^{k}T_{i, t}^{k}X_{i, t}}+\delta_{t}+\varepsilon_{i, t}`

*/

****************************************
           * PRELIMINARIES *
****************************************

set more off
set matsize 10000

local dirInput "~/Dropbox/ChicagoRAs_Reanalysis/interpolation/data/consolidated"
local dirOutput "~/Dropbox/GCP/MORTALITY/tables/Table_2/Panel_A"
cap mkdir "`dirOutput'"

use "`dirInput'/global_mortality_admin2_year_3agegrp.dta", clear

* create indicator for countries with and without age data 
gen allage_country = 1 if country == "china" | country == "india"
replace allage_country = 0 if allage_country != 1
gen agegrp_country = 1 if allage_country == 0
replace agegrp_country = 0 if agegrp_country != 1

* temperature bin variables
local BEST_bins "BEST_bin_nInf_n17C_pop BEST_bin_n17C_n12C_pop BEST_bin_n12C_n7C_pop BEST_bin_n7C_n2C_pop BEST_bin_n2C_3C_pop BEST_bin_3C_8C_pop BEST_bin_8C_13C_pop BEST_bin_18C_23C_pop BEST_bin_23C_28C_pop BEST_bin_28C_33C_pop BEST_bin_33C_Inf_pop"
local BEST_bins_lag1 "L1_BEST_bin_nInf_n17C_pop L1_BEST_bin_n17C_n12C_pop L1_BEST_bin_n12C_n7C_pop L1_BEST_bin_n7C_n2C_pop L1_BEST_bin_n2C_3C_pop L1_BEST_bin_3C_8C_pop L1_BEST_bin_8C_13C_pop L1_BEST_bin_18C_23C_pop L1_BEST_bin_23C_28C_pop L1_BEST_bin_28C_33C_pop L1_BEST_bin_33C_Inf_pop"
* precipitation interaction
local prcp_quad "precip_pop precip_pop2"
local prcp_by_adm0 "c.precip_pop#i.cntry c.precip_pop2#i.cntry"
local prcp_by_adm1 "c.precip_pop#i.adm1 c.precip_pop2#i.adm1"

cap mkdir "`dirOutput'/ster"

*save so that age shares by group can be merged back in as separate variables for the collapsed regression
tempfile collapsed
save `collapsed', replace

* merge age shares
tempfile stackedFile
reshape long deaths population deathrate adm2popshare avg_adm2popshare, i(year adm2) j(agegrp)
merge m:1 adm2 year using `collapsed', keepusing(adm2popshare31 adm2popshare32 adm2popshare33 avg_adm2popshare31 avg_adm2popshare32 avg_adm2popshare33) nogen
* fill in age shares for countries without age data

replace avg_adm2popshare31 = adm2popshare31 if allage_country == 1
replace avg_adm2popshare32 = adm2popshare32 if allage_country == 1
replace avg_adm2popshare33 = adm2popshare33 if allage_country == 1


****************************************
            * ESTIMATION *
****************************************

*** Column (1)
* global year FE
local spec1 "global_year_FE"
reghdfe deathrate `BEST_bins' if agegrp==0, absorb(`prcp_by_adm0' i.adm2#c.avg_adm2popshare31 i.adm2#c.avg_adm2popshare32 i.adm2#c.avg_adm2popshare33 i.year) cluster(adm1)
estimates save "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec1'", replace


*** Column (2)
* country-year FE
local spec2 "country_year_FE"
reghdfe deathrate `BEST_bins' if agegrp==0, absorb(`prcp_by_adm0' i.adm2#c.avg_adm2popshare31 i.adm2#c.avg_adm2popshare32 i.adm2#c.avg_adm2popshare33 i.cntry#i.year) cluster(adm1)
estimates save "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec2'", replace

*** Column (3)
* country-year-age group FE
local spec3 "country_year_age_FE"
reghdfe deathrate `BEST_bins' if agegrp==0, absorb(`prcp_by_adm0' i.adm2#c.avg_adm2popshare31 i.adm2#c.avg_adm2popshare32 i.adm2#c.avg_adm2popshare33 i.year#i.cntry#c.avg_adm2popshare31 i.year#i.cntry#c.avg_adm2popshare32 i.year#i.cntry#c.avg_adm2popshare33) cluster(adm1)
estimates save "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec3'", replace

*** Column (4)
* country-year-age group FE with state linear trends
local spec4 "country_year_age_FE_state_trends"
reghdfe deathrate `BEST_bins' if agegrp==0, absorb(`prcp_by_adm0' adm1#c.trend i.adm2#c.avg_adm2popshare31 i.adm2#c.avg_adm2popshare32 i.adm2#c.avg_adm2popshare33 i.year#i.cntry#c.avg_adm2popshare31 i.year#i.cntry#c.avg_adm2popshare32 i.year#i.cntry#c.avg_adm2popshare33) cluster(adm1)
estimates save "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec4'", replace

*** Column (5)
* state-year-age group FE
local spec5 "state_year_age_FE"
reghdfe deathrate `BEST_bins' if agegrp==0, absorb(`prcp_by_adm0' i.adm2#c.avg_adm2popshare31 i.adm2#c.avg_adm2popshare32 i.adm2#c.avg_adm2popshare33 i.year#i.adm1#c.avg_adm2popshare31 i.year#i.adm1#c.avg_adm2popshare32 i.year#i.adm1#c.avg_adm2popshare33) cluster(adm1)
estimates save "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec5'", replace

*** Column (6)
* country-year-age group FE with lagged bins
local spec6 "country_year_age_FE_lag1"
reghdfe deathrate `BEST_bins' `BEST_bins_lag1' if agegrp==0, absorb(`prcp_by_adm0' i.adm2#c.avg_adm2popshare31 i.adm2#c.avg_adm2popshare32 i.adm2#c.avg_adm2popshare33 i.year#i.cntry#c.avg_adm2popshare31 i.year#i.cntry#c.avg_adm2popshare32 i.year#i.cntry#c.avg_adm2popshare33) cluster(adm1)
estimates save "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec6'", replace


****************************************
              * OUTPUT *
****************************************

local dirInput "~/Dropbox/ChicagoRAs_Reanalysis/interpolation/data/consolidated"
local dirOutput "~/Dropbox/GCP/MORTALITY/tables/Table_2/Panel_A"
local spec1 "global_year_FE"
local spec2 "country_year_FE"
local spec3 "country_year_age_FE"
local spec4 "country_year_age_FE_state_trends"
local spec5 "state_year_age_FE"
local spec6 "country_year_age_FE_lag1"


***** Tables *****

*** Column (1) - global year FE
estimates use "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec1'"
outreg2 using "`dirOutput'/Table2_6countries_PanelA_adm0prcp.tex", tex(frag) excel replace keep(`BEST_bins') ///
          addtext("Admin2-Age FE", "YES", "Year FE", "YES", "Country-Year FE", "NO", "Country-Year-Age FE", "NO", "State time trends", "NO", "State-Year-Age FE", "NO", "1 period lags", "NO") nocons 

*** Column (2) - country-year FE
estimates use "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec2'"
outreg2 using "`dirOutput'/Table2_6countries_PanelA_adm0prcp.tex", tex(frag) excel append keep(`BEST_bins') ///
          addtext("Admin2-Age FE", "YES", "Year FE", "NO", "Country-Year FE", "YES", "Country-Year-Age FE", "NO", "State time trends", "NO", "State-Year-Age FE", "NO", "1 period lags", "NO") nocons 

*** Column (3) - country-year-age group FE
estimates use "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec3'"
outreg2 using "`dirOutput'/Table2_6countries_PanelA_adm0prcp.tex", tex(frag) excel append keep(`BEST_bins') ///
          addtext("Admin2-Age FE", "YES", "Year FE", "NO", "Country-Year FE", "NO", "Country-Year-Age FE", "YES", "State time trends", "NO", "State-Year-Age FE", "NO", "1 period lags", "NO") nocons 

*** Column (4) - country-year-age group FE with state linear trends
local spec "country_year_age_FE_state_trends"
estimates use "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec4'"
outreg2 using "`dirOutput'/Table2_6countries_PanelA_adm0prcp.tex", tex(frag) excel append keep(`BEST_bins') ///
          addtext("Admin2-Age FE", "YES", "Year FE", "NO", "Country-Year FE", "NO", "Country-Year-Age FE", "YES", "State time trends", "YES", "State-Year-Age FE", "NO", "1 period lags", "NO") nocons 

*** Column (5) - state-year-age group FE
estimates use "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec5'"
outreg2 using "`dirOutput'/Table2_6countries_PanelA_adm0prcp.tex", tex(frag) excel append keep(`BEST_bins') ///
          addtext("Admin2-Age FE", "YES", "Year FE", "NO", "Country-Year FE", "NO", "Country-Year-Age FE", "NO", "State time trends", "NO", "State-Year-Age FE", "YES", "1 period lags", "NO") nocons 

*/** Column (6) - country-year-age group FE with lagged bins
estimates use "`dirOutput'/ster/global_6countries_response_adm0prcp_`spec6'"
outreg2 using "`dirOutput'/Table2_6countries_PanelA_adm0prcp.tex", tex(frag) excel append keep(`BEST_bins') ///
          addtext("Admin2-Age FE", "YES", "Year FE", "NO", "Country-Year FE", "NO", "Country-Year-Age FE", "YES", "State time trends", "NO", "State-Year-Age FE", "NO", "1 period lags", "YES") nocons 

************ ATTENTION ************
* REPLACE THE VALUES IN COLUMN 6 WITH THOSE DISPLAYED IN THE MATRIX BELOW *
mat lags = J(12,2,0)
local bb = 0
local BEST_bins "BEST_bin_nInf_n17C_pop BEST_bin_n17C_n12C_pop BEST_bin_n12C_n7C_pop BEST_bin_n7C_n2C_pop BEST_bin_n2C_3C_pop BEST_bin_3C_8C_pop BEST_bin_8C_13C_pop BEST_bin_18C_23C_pop BEST_bin_23C_28C_pop BEST_bin_28C_33C_pop BEST_bin_33C_Inf_pop"
foreach bin in `BEST_bins' {
    local bb = `bb' + 1
    if `bb' == 8 {
    local bb = 9
    }
    lincom `bin' + L1_`bin'
    mat lags[`bb',1] = r(estimate)
    mat lags[`bb',2] = r(se)
}
matname lags estimate se, columns(1..2) explicit
matrix list lags
************ ATTENTION ************