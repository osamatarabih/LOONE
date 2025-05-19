# LOONE Data Description

## Stage data
- Water Level data download R Script: download stage data at different gauges in Lake Okeechobee and merge them in one csv file. In addition, water levels at the different WCAs locations are downloaded.
- Then using the LOONE_DATA_PREP Python Script, it uses the Stg_Sto_Ar Python script to calculate the surface area and storage associated with the water depth.
- LOONE_DATA_PREP Script put the WCA water level data in format for LOONE.

## Flows
- Using Inflow Flow data download R Script: We download flow data at all inlet locations (S191_S, S65E_S, S65EX1_S, S84_S, S154_C, S71_S, S72_S, FISHP, S308.DS, L8.441, S133_P, S127_C, S127_P, S129_C, S135_C, S2_P, S3_P, S4_P, S351_S, S352_S, S354_S, S129 PMP_P, S135 PMP_P).
- Using Outflows Flow data download R Script: we download flow data at all distributaries of Lake Okeechobee (S308.DS, S77_S, L8.441, S127_C, S129_C, S135_C, S351_S, S352_S, S354_S, INDUST, S79, S80, S2_NNR, S3, S48_S, S49_S).
- Using the LOONE_DATA_PREP Python Script, we calculate Inflows, Outflows, and Netflows.
- LO_Inflows_BK_LORS20082023, INDUST_Outflow_20082023, Netflows_acft_LORS20082023, Outflows_consd_20082023, C43RO_LORS20082023, C43RO_Monthly_LORS20082023, C44RO_LORS20082023, C44RO_Monthly_LORS20082023, SLTRIB_Monthly_LORS20082023, TotalQWCA_Obs_LORS20082023, EAA_MIA_RUNOFF_Inputs_LORS20082023 are outputs of LOONE_DATA_PREP Script using flow data downloads from DBHYDRO.

## Weather Data
- Lake Weather data download R Script is used to download all meteorological data (RF, ET, Wind, RADT, RADP, Water Temp, Air Temp) at different stations in the lake.
- Rainfall data were downloaded and merged into one file and average daily rainfall over the lake was determined using the R Script.
- ET data were downloaded and merged into one file and average daily rainfall over the lake was determined using the R Script.
- Using the LOSMB_Updated Python Script: The Storage deviation is calculated using Rainfall, ET, Netflows, Outflows, S77_In, and S308_In data calculated using the LOONE_DATA_PREP Script.
- For Water Temperature, I manually merge the data at all stations into one file and calculate the average of all values each day.
- We used regression relations developed by our team to predict water temperature function of air temperature because water temperature measurements stopped around 2012 at most stations and around 2020 at one station only (LZ40).
- Then, I used Python functions to fill all the gaps in the Water Temperature data with predicted water temperature (i.e., function of air temperature).
- The Water Temperature data were used to calculate kinematic viscosity (nu) using the Kinematic_Visc Python Script.
- Wind Speed data were averaged using the LOONE_DATA_PREP Script.
- Then, using the Wind_Induced_Waves Python Script, the shear stress is calculated.
- The output of this Script is File (WindShearStress_LORS20082023).

## Water Quality Data
- Using Inflow Water Quality data download R Script: We download PHOSPHATE, TOTAL AS P, AMMONIA-N, NITRATE+NITRITE-N, CHLOROPHYLL-A, CORRECTED, CHLOROPHYLL-A(LC) for main inflow stations (S191', 'S65E', 'S84', 'S154', 'S71', 'S72', 'S4', 'FECSR78', 'S308C', 'CULV10A', 'S133', 'S127', 'S135').
- Data Interpolation
- Using Data_Interpolation Python Script, we interpolate all the water quality in the previous step to daily values. Note: Needs to run all stations at the same time and interpolate all simultaneously.
- I did interpolate NITRATE+NITRITE-N, AMMONIA-N, PHOSPHATE, TOTAL AS P, PHOSPHATE, ORTHO AS P, CHLOROPHYLL-A, CORRECTED, CHLOROPHYLL-A(LC) for all stations inside the lake and for all inflow stations.
- Using LOONE_DATA_PREP Python Script average daily of TP, OP, NO, NH4, Chla are calculated in Lake Okeechobee including all available stations.
- Using LOONE_DATA_PREP Python Script, we calculate TP loads at all inflow stations, Then Total External Loads and Total External Loads 3MLag are calculated.

## NOTES
- Files CE_SLE_turns_inputs_LORS20082023, Estuary_needs_water_Input_LORS20082023, Water_dmd_LORS20082023, as well as columns LOSAdmd and LOSAsup in the SFWMM_Daily_Outputs_LORS20082023 were filled using replicated data from the past!
- PhotoPeriod Python Script is used to determine daily photoperiod for one year in Lake Okeechobee and then this one year data is replicated to fill in all years of simulation.
- Seasonal_LONINO_LORS20082023 and Multi_Seasonal_LONINO_LORS20082023 were developed manually by downloading monthly USACE Reports and fill missing data.
