# LOONE
The Lake Operation Optimization of Nutrient Exports (LOONE) is a comprehensive water balance-nutrient-optimization model that comprises three coupled modules: A water balance module that simulates the water balance and operations of a reservoir, a nutrient module that simulates nutrient (phosphorus) dynamics in the water column, and an optimization engine that optimizes a reservoir’s releases into its distributaries with the objective of nutrient export minimization and/or water deficit minimization. Python 3 was chosen to develop the code because of its many high-quality libraries and a powerful community support. 

## Installation
```bash
pip install loone
```

### Development Installation
```bash
git clone <this repository>
cd ./LOONE
pip install -e .
```
 
## How to Run LOONE?
```python
"""
Data prep
    1. Add all required data to the workspace directory.
    2. Add a config.yaml file following the example below with the correct variables and file names for required data.
"""
from loone.loone_q import LOONE_Q
from loone.loone_nut import LOONE_NUT
from loone.loone_wq import LOONE_WQ


LOONE_Q(
    workspace="/path/to/workspace",
    p1=0,
    p2=0,
    s77_dv=0,
    s308_dv=0,
    tp_lake_s=0,
)

LOONE_NUT(
    workspace="/path/to/workspace",
    out_file_name="loone_nut_outputs.csv",
    loads_external_filename="lo_external_loads.csv",
    flow_df_filename="flow_df.csv",
    forecast_mode=True,
)

LOONE_WQ(workspace="/path/to/workspace")
```


### Example configuration file
```yaml
# LOONE Configuration

# predefined variables
schedule: "LORS20082023"
sim_type: 0  # 0:Scenario_Simulation 1:Optimization_Validation 2:Optimization 3:Loone_Scenarios_App_Simulation
start_year: 2008
end_year: 2023
start_date_entry: [2008, 1, 1]
beg_date_cs_entry: [2008, 1, 1]
end_date_entry: [2023, 3, 31]
end_date_tc: [2023, 4, 1]
month_n: 183
opt_new_tree: 1  # if New Tree Decision is used enter 1 else enter 0.
code: 6
multiplier: 100
tci: 1
opt_net_inflow: 2
net_inf_const: 0
start_stage: 10.268
beg_stage_cs: 10.268
opt_los_admd: 1
mult_losa: 100
opt_losa_ws: 1  # the option for LOSA Daily Supply where 1: Calculated Function of WSM, 2: Set values to zeros.
opt_dec_tree: 1  # if Tree Decision is used enter 1 else enter 0.
zone_c_met_fcast_indicator: 1  # 0 to use the same tree classifications as SLONINO or 1 to use the same tree classifications as LORS2008 and SFWMM.
wca3a_reg_zone: "ERTP:TopE"
wca3a_offset: 0
wca328_min: 7.5
wca3_nw_min: 11
wca217_min: 11.1
opt_wca_limit_wsa: 2
cs_flag: 1
pls_day_switch: 0  # 0: pulse day counter continues to 10 even if release level increases, 1: pulse day counter is set to zero if release level increases during the 10-day pulse.
max_qstg_trigger: 20  # the maximum stage trigger for maximum discharge if Trib_cond. = XWet.
opt_qreg_mult: 0  # option for using multipliers 0: don't use, 1: apply only during dry season (Nov-May), 2: apply year-round.
alternate_high_qyrs: 0
option_s80_baseflow: 0
s308_bk_const: 1
s308_bk_thr: 14.5
opt_s308: 1
s308_rg_const: 1
option_reg_s77_s308: 0
s80_const: 1
opt_outlet1_dsrg: 0
thc_threshold: 2
low_chance: 50
opt_l_chance_line: 1
opt_date_targ_stg: 1
opt_sal_fcast: 3
ce_sal_threshold: 5
late_dry_season_option: 0
opt_no_ap_above_bf_sb: 1
opt_adap_prot: 1
opt_ceews_lowsm: 0
opt_thc_byp_late_ds: 1
apcb1: 100
apcb2: 100
apcb3: 100
apcb4: 100
cal_est_ews: 300
outlet1_usews_switch: 1
outlet1_usbk_switch: 1  # option for S77BK simulation 0: Use input data or 1: Simulate with LOONE.
outlet1_usbk_threshold: 11.1
option_s77_baseflow: 0  # 0: baseflow supplements daily C43RO, 1: baseflow supplements monthly C43RO.
outlet1_usreg_switch: 1
outlet1_ds_switch: 1
max_cap_reg_wca: 4000
multiplier_reg_wca: 1
option_reg_wca: 2
constant_reg_wca: 400
max_cap_reg_l8_c51: 500
multiplier_reg_l8_c51: 1
option_reg_l8_c51: 2
constant_reg_l8_c51: 200
et_switch: 0
opt_wsa: 0  # Options for Lake O water supply augmentation (WSA) operation (0 = no WSA,1 = use flat trigger stages to activate WSA operation,2 = trigger stages defined using offsets from LOWSM WST line to activate WSA operation)
wsa_thc: 2
wsa_trig1: 12.5
wsa_trig2: 11.5
wsa_off1: 0.5
wsa_off2: 0
mia_cap1: 860
mia_cap2: 1720
nnr_cap1: 900
nnr_cap2: 1800
option_stage: 0
# Water demand cutback for each WSM Zone
z1_cutback: 0.15
z2_cutback: 0.3
z3_cutback: 0.45
z4_cutback: 0.6

dstar_b: 99
dstar_c: 99
dstar_d3: 99
dstar_d2: 99
dstar_d1: 99
astar_b: 1
astar_c: 1
astar_d3: 1
astar_d2: 1
astar_d1: 1
bstar_s77_b: 1
bstar_s77_c: 1
bstar_s77_d3: 0.5
bstar_s77_d2: 0.5
bstar_s77_d1: 0.5
bstar_s80_b: 1
bstar_s80_c: 1
bstar_s80_d3: 0.5
bstar_s80_d2: 0.5
bstar_s80_d1: 0.5

# data
sfwmm_daily_outputs: "SFWMM_Daily_Outputs.csv"
wsms_rsbps: "WSMs_RSBPs.csv"
losa_wkly_dmd: "LOSA_wkly_dmd.csv"
trib_cond_wkly_data: "Trib_cond_wkly_data.csv"
seasonal_lonino: "Seasonal_LONINO.csv"
multi_seasonal_lonino: "Multi_Seasonal_LONINO.csv"
netflows_acft: "Netflows_acft.csv"
water_dmd: "Water_dmd.csv"
rf_vol: "RFVol.csv"
et_vol: "ETVol.csv"
c44ro: "C44RO.csv"
c43ro: "C43RO.csv"
basin_ro_inputs: "Basin_RO_inputs.csv"
c43ro_monthly: "C43RO_Monthly.csv"
c44ro_nonthly: "C44RO_Monthly.csv"
sltrib_monthly: "SLTRIB_Monthly.csv"
s77_regulatory_release_rates: "S77_RegRelRates.csv"
s80_regulatory_release_rates: "S80_RegRelRates.csv"
ce_sle_turns_inputs: "CE_SLE_turns_inputs.csv"
pulses_inputs: "Pulses_Inputs.csv"
june_1st_lake_stage_below_11ft: "Chance of June 1st Lake stage falling below 11.0ft.csv"
may_1st_lake_stage_below_11ft: "Chance of May 1st Lake stage falling below 11.0ft.csv"
estuary_needs_water_input: "Estuary_needs_water_Input.csv"
eaa_mia_ro_inputs: "EAA_MIA_RUNOFF_Inputs.csv"
storage_deviation: "Storage_Dev.csv"
calibration_parameters: "Cal_Par.csv"

# tp variables regions
z_sed: 0.05  # m
per_h2o_m: 0  # 85 #%
per_h2o_s: 0  # 20 #%
per_h2o_r: 0  # 20  #%
per_h2o_p: 0  # 85 #%
n_per: 0.43
s_per: 0.57

bulk_density_m: 0.15  # g/cm3
bulk_density_s: 1.213  # g/cm3
bulk_density_r: 1.213  # g/cm3
bulk_density_p: 0.14  # g/cm3

particle_density_m: 1.2  # g/cm3
particle_density_s: 2.56  # g/cm3
particle_density_r: 2.56  # g/cm3
particle_density_p: 1.2  # g/cm3

a_mud_n: 377415128  # m2 in 1988!
a_mud_s: 394290227  # m2 in 1988!
a_sand_n: 237504380  # m2 in 1988!
a_sand_s: 117504905  # m2 in 1988!
a_rock_n: 17760274  # m2 in 1988!
a_rock_s: 141327951  # m2 in 1988!
a_peat_n: 97497728  # m2 in 1988!
a_peat_s: 301740272  # m2 in 1988!

v_burial_m: 0.0000003  # 1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)
v_burial_s: 0.0000003  # 1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)
v_burial_r: 0.0000003  # 1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)
v_burial_p: 0.0000003  # 1.0e-05 #(m/day)#0.00017333#(m/month)# 0.00208 (m/yr)

nondominated_sol_var: "nondominated_Sol_var.csv"

wca_stages_inputs: "WCA_Stages_Inputs.csv"
lo_inflows_bk: "LO_Inflows_BK.csv"
sto_stage: "Average_LO_Storage_3MLag.csv"
wind_shear_stress: "WindShearStress.csv"
nu: "nu.csv"

outflows_observed: "Flow_df_3MLag.csv"
```

## Case Study:
LOONE was used to simulate operations and phosphorus mass balance of Lake Okeechobee, the largest reservoir by surface area in the US. In addition to its dimensions, we chose Lake Okeechobee as a case study for multiple reasons. First, it is a multi-inlet and multi-outlet reservoir with a complex system of pumps and locks operated with a regulation schedule, which allowed us to evaluate performance of the model in a complex hydrologic system. Second, it is a nutrient impaired lake known to export large amount of nutrients to its surrounding water bodies (Tarabih and Arias, 2021; Walker, 2000), thus LOONE could aid evaluating impacts of lake regulations on the nutrient status of the regional system.

## Data Description:
See : [Data Description](loone/data/data_description.md)

## Data Requirements:
LOONE was used to simulate Lake Okeechobee regulatory releases into St. Lucie Canal, and Caloosahatchee River, meanwhile prescribed flows were used for West Palm Beach Canal, North New River Canal/Hillsboro Canal, Miami Canal, and L-8 Canal as well as water supply flows utilizing the continuity equation and Lake Okeechobee rule curves for the study period. We simulated three different Lake Okeechobee schedules during the study period (1991-2018): RUN25 (1991-1999), WSE (2000-2007), and 2008 LORS (2008-2018). 
LOONE was used to design optimal releases of Lake Okeechobee into the Caloosahatchee River and St. Lucie Canal, with the goal of demonstrating an operational schedule that can minimize pollutant exports into the estuaries while minimizing LOSA water deficits. 

| Data type              | Explanation                                                | Time step | File name                      | Data source            |
|------------------------|------------------------------------------------------------|-----------|--------------------------------|------------------------|
| Tributary Condition    | Net Rainfall, Tributary Flow, Palmer Index, Net inflows    | Weekly    | Trib_cond_wkly_data_xxx       | Rainfall, Tributary flow, and net inflows (DBHYDRO) |
| Palmer Index           | NOAA                                                       | Monthly   | Palmer_Index_xxx              | USACE Monthly Reports  |
| Seasonal LONINO        | Seasonal Lake Okeechobee Net Inflow Outlooks               | Monthly   | Seasonal_LONINO_xxx           | USACE Monthly Reports  |
| Multi Seasonal LONINO  | Multi Seasonal Lake Okeechobee Net Inflow Outlooks         | Monthly   | Multi_Seasonal_LONINO_xxx     | USACE Monthly Reports  |
| Net Inflows            | Net Inflows = All tributary inflows – non simulated outflows | Daily     | NetFlows_acft_xxx             | DBHYDRO                |
| Water demand           | LOSA water demand                                          | Daily     | Water_dmd_xxx                 | SFWMD Reports          |
| Rainfall               | Rainfall Volume                                            | Daily     | RF_Volume_xxx                 | DBHYDRO                |
| Evapotranspiration     | ET Volume                                                  | Daily     | ETVol_xxx                     | DBHYDRO                |
| C44 Runoff             | St Lucie Watershed Runoff                                  | Daily     | C44RO_xxx                     | DBHYDRO                |
| C43 Runoff             | Caloosahatchee Watershed Runoff                            | Daily     | C43RO_xxx                     | DBHYDRO                |
| EAA_MIA_Runoff         | Daily flow data for Miami Canal at S3, NNR at S2_NNR, WPB at S352, S2 pump, and S3 pump. | Daily | EAA_MIA_RUNOFF_Inputs_xxx | DBHYDRO |
| Storage Deviation      | Storage deviation between simulated storage using observed outflows and observed storage to account for unreported outflows. | Daily | Storage_Dev_xxx | DBHYDRO |
| External Loads         | Phosphorus loads into Lake Okeechobee from the tributaries | Daily     | LO_External_Loadings_3MLag_xxx | DBHYDRO |
| Lake Inflows           | Lake Okeechobee inflows from all the tributaries as well as back flows | Daily | LO_Inflows_BK_xxx | DBHYDRO |
| Wind Shear Stress      | Wind shear stress function of wind speed                    | Daily     | WindShearStress_xxx           | Calculated             |
| Wind Speed             | Mean wind speed                                            | Daily     | Mean_WindSpeed_xxx            | DBHYDRO                |
| Kinematic viscosity    | Kinematic viscosity of Lake Okeechobee water column function of Water Temperature | Daily | nu_xxx | DBHYDRO |
| Water Temperature      | Water column Temperature                                   | Daily     | LZ40_T_xxx                    | DBHYDRO                |

