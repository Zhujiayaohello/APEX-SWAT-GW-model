APEX-SWAT-GW model with R quick start manual
================


Introduction
------------
The linked APEX-SWAT-GW model was constructed to transfer variable values between the APEX and SWAT-GW simulation subroutines.The linking of these two models is forced by transferring the PRK (APEX variable) to Wrchrg (SWAT-GW variable).


Steps for linking APEX and SWAT-GW
-----------------------
``` r
# 1)library(data.table)
   library(tidyverse)

# 2)set the path
path <- 'E:/APEX-SWAT-GW/APEX-SWAT-GW with R'
setwd(path)
getwd()

# 3)read ****.MSA file in APEX folder
sit1 <-
  read.table(
    'SITE1.MSA',
    skip = 10,
    comment.char = "",
    header = F,
    fill = T
  )
sit1r <-
  rename(sit1, c(
    "GSub" = "V2",
    "Year" = "V3",
    "class" = "V5",
    "value" = "V18"
  ))
sit1r <- sit1r[, c("GSub", "Year", "class", "value")]

# 4)filter required parameters
sit1_wid <-
  pivot_wider(sit1r, names_from = class, values_from = value)
sit1_s <- select(sit1_wid, GSub, Year, PRCP, IRGA, ET, WYLD)
fwrite(sit1_s, "sit1.csv")

# 5)scientific notation is not allowed
op <- options(scipen = 99999)
sit_start <- fread("sit1.csv")

# 6)calculate PRK
sit1_PRK <- mutate(sit_start, PRK = PRCP + IRGA - ET - WYLD)
sit1_PRK <- select(sit1_PRK, GSub, Year, PRK)
head(sit1_PRK)

# 7)filter required subbasins
#sit1_PRK <- filter(sit1_PRK, GSub == 3 | GSub == 4| GSub == 5| GSub == 6| GSub == 8| GSub == 9| GSub == 10| GSub == 13| GSub == 15| GSub == 17| GSub == 18| GSub == 21)
#sit1_PRK<- filter(sit1_PRK, Year != 1990 & Year != 1991)

# 8)read **gw.CSV from SWAT-GW folder
sit1_gw <- read.table('zygw.csv', header = T, sep = ',')
zy_merge <- merge(sit1_gw, sit1_PRK,
                  by = intersect(names(sit1_gw), names(sit1_PRK)))
head(zy_merge)

# 9)calculate groundwater table variation (H)
zy_S <-
  mutate(zy_merge, S = (PRK + LARCHRG - DA_RCHG - SA_IRR - Qgw - REVAP))
zy_H <- mutate(zy_S, H = (S / u) / 1000)
head(zy_H)

# 10)output calculation results **_H.CSV
fwrite(zy_H, "zy_H.csv")
```

APENDIX
-------

This section provides the name of different parameters and their description.

| Parameters                | Description                                                                                                                         |
|:--------------------------|:------------------------------------------------------------------------------------------------------------------------------------|
| GSub                      | Subbasin number as it appears inside SWAT-GW **gw.CSV file.                                                                         |
| DA_RCHG                   | Wdeep:Amount of water entering deep aquifer from root zone as it appears inside SWAT-GW **gw.CSV file (mm).                         |
| SA_IRR                    | Wpump:Amount of water removed from shallow aquifer for irrigation as it appears inside SWAT-GW **gw.CSV file (mm).                  |
| REVAP                     | Wrevap:Water in shallow aquifer returning to root zone as it appears inside SWAT-GW **gw.CSV file (mm).                             |
| LARCHRG                   | Wlarchrg:Amount of lateral recharge from mountainous area entering into the aquifer as it appears inside SWAT-GW **gw.CSV file (mm).|
| Qgw                       | Groundwater discharge into reach as it appears inside SWAT-GW **gw.CSV file (mm).                                                   |
| u                         | specific yield as it appears inside SWAT-GW **gw.CSV file (m3 m-3).                                                                 |
| PRK                       | Percolation as it calculate inside R (mm).                                                                                          |
| PRCP                      | precipitation as it appears inside APEX ****.MSA file (mm).                                                                         |
| IRGA                      | Irrigation as it appears inside APEX ****.MSA file (mm).                                                                            |
| ET                        | evapotranspiration as it appears inside APEX ****.MSA file (mm).                                                                    |
| WYLD                      | water yield as it appears inside APEX ****.MSA file (mm).                                                                           |
| S                         | the amount of water stored in the shallow aquiferas it appears inside Output **_H.CSV file (mm).                                    |
| H                         | groundwater table variation as it appears inside Output **_H.CSV file (m).                                                          |



====================================================================================================================================================================


APEXSENSUN quick start manual
================


Introduction
------------

APEXSENSUN is a package in R for performing uncertainty and sensitivity analysis (SA) for the APEX model.The improved APEXSENSUN replace the read-in DWS file with a read-in SAD file, which can print the daily ET at each subarea, and made it more flexible in performing sensitivity with the function of selecting different subareas and different crop covers.


Example folder
--------------

An example folder containing an APEX project and other inputs is available for users to test the package. The rest of this manual provides details of implementing an SA project using the accompanying example folder which can be created through a call to:

``` r
# Loading APEXSENSUN package in R:
  library(APEXSENSUN)

# Opening Modified APEXSENSUN R source code and using it
  
# Creating a copy of tutorial folder inside the working directory
  getExampleFolder()
```

Steps for performing SA
-----------------------

After loading APEXSENSUN and generating a copy of the example folder, the following four steps, described in the next sections should be followed for performing SA.

``` r
# 1) Generating a list object with a predefined structure compatible to APEXSENSUN:
     globalInput <- APEXSENSUN::inputGen()

# 2) Setting the required inputs (e.g. uncertainty boubds, SA method, sample size, ...)
  #
  # Setting uncertainty bounds:
    globalInput$apexPARM$Root_growth_soil[1] = 0.15
    globalInput$apexPARM$Root_growth_soil[2] = 0.2
  
    globalInput$apexPARM$Soil_water_limit[1] = 0
    globalInput$apexPARM$Soil_water_limit[2] = 1
  
    globalInput$apexPARM$Soil_evap_coeff[1] = 1.5
    globalInput$apexPARM$Soil_evap_coeff[2] = 2.5
  
    globalInput$apexPARM$Soil_evap_plant_cover[1] = 0
    globalInput$apexPARM$Soil_evap_plant_cover[2] = 0.5
  
    globalInput$apexPARM$Runoff_CN_int_abs[1] = 0.05
    globalInput$apexPARM$Runoff_CN_int_abs[2] = 0.4
  
    globalInput$apexPARM$Max_rain_intercept[1] = 0
    globalInput$apexPARM$Max_rain_intercept[2] = 15

  
# SA method and sample size:
    globalInput$gsaType <- "SRRC"
    globalInput$sampleSize <- 1000

# 3) Performing Monte Carlo simulation using the setting in globalInput:
     input4SA <- mc4APEX(globalInput)
    
# 4) Calculation of sensitivity indices:
     sa4APEX(globalInput,input4SA = input4SA)
```


APENDIX
-------

This section provides 3 main tables containing the name of different parameters and their description.

Table: Table 1. General inputs

| Parameters                | Description                                                                                                                         |
|:--------------------------|:------------------------------------------------------------------------------------------------------------------------------------|
| sampleSize                | Sample size or length of the discretization of parameter space.                                                                     |
| captionVarSim             | Simulated variable name as it appears inside .DWS APEX file.                                                                        |
| captionVarObs             | Observed variable names as appears inside observed file.                                                                            |
| startDate                 | Start date for analysis with format: YYYY MM DD e.g., 2002 01 01                                                                    |
| endDate                   | End date for analysis with format: YYYY MM DD e.g., 2003 01 01                                                                      |
| labelAPEXExe              | APEX executable file name excluding file's extension.                                                                               |
| labelWatershedParam       | APEX PARM file name excluding file's extension.                                                                                     |
| labelControlParam         | APEXCONT file name excluding file's extension.                                                                                      |
| labelOutputVariableAWP    | APEX .AWP file name excluding file's extension.                                                                                     |
| labelOutputVariableACY    | APEX .ACY file name excluding file's extension.                                                                                     |
| labelOutputVariableDWS    | APEX .DWS file name excluding file's extension.                                                                                     |
| labelObservedVar.txt      | Observed file name containing observed time series.                                                                                 |
| backUpPARM0806.dat        | Path to original APEX file containing PARM parameters.                                                                              |
| backUpAPEXCONT.dat        | Path to original APEX file containing APEXCONT parameters.                                                                          |
| folderPathProject         | Path to folder containing APEX model.                                                                                               |
| folderPathRCodes          | No need to set!                                                                                                                     |
| folderPathObserved        | Path to folder containing observed data file.                                                                                       |
| folderPathGsaOutputs      | Path to folder storing SA results.                                                                                                  |
| storeFolderPathWatershed  | Path to folder storing generated PARM files for Monte Carlo runs.                                                                   |
| storeFolderPathControl    | Path to folder storing generated APEXCONT file for Monte Carlo runs.                                                                |
| calculatedOutputFolderAWP | Path to folder storing calculated .AWP files for Monte Carlo runs.                                                                  |
| calculatedOutputFolderACY | Path to folder storing calculated .ACY files for Monte Carlo runs.                                                                  |
| calculatedOutputFolderDWS | Path to folder storing calculated .DWS files for Monte Carlo runs.                                                                  |
| gsaType                   | Type of SA method: MORRIS, SRC, SRRC, SOBOL, SOBOL2002, SOBOL2007, SOBOLEFF, SOBOLJANSEN, SOBOLMARA, SOBOLMARTINEZ, FAST99, KSTEST. |

Table: Table 2. SA-specific parameters                                                                                        

| Parameters       | Description                                                                                                  |
|:-----------------|:-------------------------------------------------------------------------------------------------------------|
| morrisRFactor    | an integer representing design repetition number (i.e. the number of elementary effect computed per factor). |
| morrisLevels     | an integer specifying the number of levels of the design in OAT (Once At a Time) design.                     |
| sobolOrder       | an integer representing maximum order in the ANOVA decomposition in Sobol method.                            |
| ksTestPerf       | A performance function type for KSTEST method. Available options are: NASH, RMSE, PBIAS                      |
| ksTestThreshold  | Threshold value for performance function for determining behavioral from non-behavioral simulations.         |
| ksTestSigmaLevel | Significance level used in KSTEST                                                                            |

Table: Table 3. APEX model parameters located inside PARM****.dat file                                                                                

| Parameters                            | Description                                                                         |
|:--------------------------------------|:------------------------------------------------------------------------------------|
| Crop\_canopy\_PET                     | Crop canopy-PET                                                                     |
| Root\_growth\_soil                    | Root growth-soil strength                                                           |
| Water\_stress\_harvest                | Water stress-harvest index                                                          |
| Water\_storage\_N                     | Water storage N leaching                                                            |
| Soil\_water\_limit                    | Soil water lower limit                                                              |
| Winter\_dormancy                      | Winter dormancy                                                                     |
| N\_fixation                           | N fixation                                                                          |
| Soluble\_P\_runoff                    | Soluble phosphorus runoff coefficient                                               |
| Pest\_damage\_moisture                | Pest damage moisture threshold                                                      |
| Pest\_damage\_cover                   | Pest damage cover threshold                                                         |
| Moisture\_req\_seed\_germ             | Moisture required for seed germination                                              |
| Soil\_evap\_coeff                     | Soil evaporation coefficient                                                        |
| Wind\_erod\_coeff                     | Wind erodibility coefficient                                                        |
| Nitrate\_leac\_ratio                  | Nitrate leaching ratio                                                              |
| Runoff\_CN\_Adj\_parm                 | Runoff CN Residue Adjustment Parameter                                              |
| Expand\_CN\_ret\_parm                 | Expands CN retention parameter                                                      |
| Soil\_evap\_plant\_cover              | Soil evaporation – plant cover factor                                               |
| Sedim\_rout\_exponent                 | Sediment routing exponent                                                           |
| Sedim\_rout\_coeff                    | Sediment routing coefficient                                                        |
| Runoff\_CN\_int\_abs                  | Runoff curve number initial abstraction                                             |
| Soluble\_C\_adsorp\_coeff             | Soluble Carbon adsorption Coefficient                                               |
| CN\_retention\_frozen\_soil           | Reduces NRCS Runoff CN Retention Parameter for Frozen Soil                          |
| Harg\_equation\_parm                  | Hargreaves PET equation coefficient                                                 |
| Pest\_leach\_ratio                    | Pesticide leaching ratio                                                            |
| Expo\_coeff\_rainfall                 | Exponential coefficient used to account for rainfall intensity on curve number      |
| Matur\_frac\_spring                   | Fraction of maturity at spring growth initiation                                    |
| CEC\_effect\_nitrification            | CEC effect on nitrification & volatilization                                        |
| N\_fixation\_limit                    | Upper Nitrogen Fixation limit                                                       |
| Biological\_mix\_efficiency           | Biological mixing efficiency                                                        |
| Soluble\_P\_exponent                  | Soluble phosphorus runoff exponent                                                  |
| Max\_depth\_bio\_mixing               | Maximum depth for biological mixing                                                 |
| OrgP\_loss\_exponent                  | Organic P loss exponent                                                             |
| MUST\_coeff                           | Coefficient in MUST EQ                                                              |
| Harg\_PET\_exponent                   | Hargreaves PET equation exponent                                                    |
| Denit\_soil\_threshold                | Denitrification soil-water threshold                                                |
| Daily\_denit\_limit                   | Upper Limit of Daily Denitrification rate                                           |
| SWAT\_delivery\_ratio\_exponent       | Exponent in Delivery Ratio for SWAT Output                                          |
| Water\_stress\_coeff                  | Water stress weighting coefficient                                                  |
| Puddling\_sat\_conduct                | Puddling Saturated conductivity                                                     |
| Groundwater\_stor\_threshold          | Groundwater storage threshold                                                       |
| Root\_temp\_stress\_exponent          | Plant root temperature stress exponent                                              |
| SCS\_index\_coeff                     | SCS curve number index coefficient                                                  |
| Plow\_depth                           | Plow layer depth                                                                    |
| CN\_retention\_param                  | Upper Limit of Curve Number Retention Parameter                                     |
| sediment\_rout\_travel\_coeff         | Sediment routing travel time coefficient                                            |
| RUSLE\_c\_factor\_res                 | RUSLE C-factor coefficient                                                          |
| RUSLE\_c\_factor\_height              | RUSLE C-factor coefficient                                                          |
| Climate\_stress\_factor               | Adjusts climatic stress factor                                                      |
| Max\_rain\_intercept                  | Maximum rainfall interception by plant canopy                                       |
| Rain\_intercept\_coeff                | Rainfall interception coefficient                                                   |
| Water\_stor\_residue\_coeff           | Water stored in litter (residue) coefficient                                        |
| Tillage\_residue\_decay\_rate\_coeff  | Exponential coefficient in EQUATION expressing tillage effect on residue decay rate |
| Microbial\_soil\_depth\_coeff         | Coefficient in oxygen EQUATION used in modifying microbial activity with soil depth |
| N\_enrich\_coeff                      | N enrichment ratio coefficient for routing                                          |
| N\_enrich\_rout\_exponent             | N enrichment ratio exponent for routing                                             |
| Fraction\_destroyed\_burn             | Fraction destroyed by burn operation                                                |
| P\_enrich\_rout\_coeff                | P enrichment ratio coefficient for routing                                          |
| P\_enrich\_rout\_exponent             | P enrichment ratio exponent for routing                                             |
| P\_move\_evap\_coeff                  | P upward movement by evaporation coefficient                                        |
| Max\_days\_grazed\_rotation           | Maximum number of days a pasture is grazed before rotation                          |
| Soil\_water\_up\_flow\_limit          | Soil water Upward Flow Limit                                                        |
| Manure\_erosion\_equation\_coeff      | Manure erosion equation coefficient                                                 |
| N\_enrich\_ratio\_delivery            | N Enrichment Ratio for Delivery to SWAT                                             |
| Dust\_distribution\_coeff             | Dust distribution coefficient                                                       |
| RUSLE2\_trans\_capacity               | RUSLE2 transport capacity parameter                                                 |
| RUSLE2\_trans\_capacity\_threshold    | RUSLE2 threshold transport capacity coefficient                                     |
| Dust\_distribution\_exponent          | Dust distribution dispersion exponent                                               |
| Manure\_erosion\_exponent             | Manure erosion exponent                                                             |
| Microbial\_top\_soil\_coeff           | Coefficient adjusts microbial activity function in the top soil layer               |
| Microbial\_decay\_coeff               | Microbial decay rate coefficient                                                    |
| Manure\_erosion\_coeff                | Manure erosion coefficient                                                          |
| Volt\_nitrification\_partition\_coeff | Volatilization/nitrification partitioning coefficient                               |
| Hydrograph\_dev\_param                | Hydrograph development parameter                                                    |
| Partition\_N\_flow\_groundwater       | Partitions Nitrogen flow from groundwater                                           |
| P\_enrich\_ratio\_deliver\_SWAT       | P Enrichment Ratio for Delivery to SWAT                                             |
| Stand\_dead\_fall\_rate\_coeff        | Standing Dead fall rate coefficient                                                 |
| Runoff\_2\_delay\_pest                | Runoff amount to delay pest application                                             |
| Soil\_water\_2\_delay\_tillage        | Soil water value to delay tillage                                                   |
| Auto\_mov\_lower\_limit               | Auto mow lower limit                                                                |
| Nitrification\_vol\_upper\_limit      | Upper Limit of Nitrification-Volatilization                                         |
| Tech\_coeff                           | Technology Coefficient                                                              |
| Drainage\_lateral\_conduct            | Estimates drainage system lateral hydraulic conductivity                            |
| P\_flux\_labile\_active\_coeff        | Coefficient regulating P flux between labile and active pool                        |
| P\_flux\_active\_stable\_coeff        | Coefficient regulating P flux between active and stable pool                        |
| N\_salt\_evap\_coeff                  | Nitrogen and Salt Upward movement by evaporation coefficient                        |
| Water\_table\_recession\_coeff        | Water table recession coefficient                                                   |
| Water\_table\_move\_limit             | Limits daily water table movement                                                   |
| Water\_table\_recession\_exponent     | Water table recession                                                               |
| Subsurface\_flow\_factor              | Subsurface flow factor                                                              |
| Flood\_evap\_limit                    | Flood Evaporation Limit                                                             |
| Runoff\_adj\_link                     | Runoff Volume Adjustment for Direct Link                                            |
| Water\_erosion\_threshold             | Water Erosion Threshold                                                             |
| Wind\_erosion\_threshold              | Wind Erosion Threshold                                                              |
| Crop\_stress\_temp\_exponent          | Exponent of Crop Stress Temperature function                                        |
| Soluble\_P\_leach\_KD                 | Soluble Phosphorus Leaching KD value                                                |
| Unknown1                              | ---                                                                                 |
| Unknown2                              | ---                                                                                 |
| Unknown3                              | ---                                                                                 |
| Irrigation\_cost                      | Cost of Irrigation Water                                                            |
| Lime\_cost                            | Cost of Lime                                                                        |
| Fuel\_cost                            | Cost of Fuel                                                                        |
| Labor\_cost                           | Cost of Labor                                                                       |
| Unknown4                              | ---                                                                                 |

