# LOONE
Lake Operation Optimization of Nutrient Exports
LOONE is a comprehensive water balance-nutrient-optimization model that comprises three coupled modules: A water balance module that simulates the water balance and operations of a reservoir, a nutrient module that simulates nutrient (phosphorus) dynamics in the water column, and an optimization engine that optimizes a reservoir’s releases into its distributaries with the objective of nutrient export minimization and/or water deficit minimization. Python 3 was chosen to develop the code because of its many high-quality libraries and a powerful community support. 
#Data Requirements:
1- Inflow time series 
2- Phosphorus mass load In
3- Rainfall time series
4- Evapotranspiration timer series
5- Sediment properties of the bottom sediment layer
6- Elevation-Storage and Elevation-Surface area for the reservoir
7- Reservoir operation guide

#How to Run LOONE?
1- Make sure all required data are provided in the data folder.
2- Configure the model main parameters via Model_Config.py where simulation type is identified and working directory path is identified.
3- Main LOONE variables are identified (e.g., start and end date of simulation, initial elevation of the lake, options to either simulate specific regulatory discharges or to use measured flow values, maximum capacity of different outlet control structures, etc.).

#Case Study:
LOONE was used to simulate operations and phosphorus mass balance of Lake Okeechobee, the largest reservoir by surface area in the US. In addition to its dimensions, we chose Lake Okeechobee as a case study for multiple reasons. First, it is a multi-inlet and multi-outlet reservoir with a complex system of pumps and locks operated with a regulation schedule, which allowed us to evaluate performance of the model in a complex hydrologic system. Second, it is a nutrient impaired lake known to export large amount of nutrients to its surrounding water bodies (Tarabih and Arias, 2021; Walker, 2000), thus LOONE could aid evaluating impacts of lake regulations on the nutrient status of the regional system.
LOONE was used to simulate Lake Okeechobee regulatory releases into St. Lucie Canal, and Caloosahatchee River, meanwhile prescribed flows were used for West Palm Beach Canal, North New River Canal/Hillsboro Canal, Miami Canal, and L-8 Canal as well as water supply flows utilizing the continuity equation and Lake Okeechobee rule curves for the study period. We simulated three different Lake Okeechobee schedules during the study period (1991-2018): RUN25 (1991-1999), WSE (2000-2007), and 2008 LORS (2008-2018). 
LOONE was used to design optimal releases of Lake Okeechobee into the Caloosahatchee River and St. Lucie Canal, with the goal of demonstrating an operational schedule that can minimize pollutant exports into the estuaries while minimizing LOSA water deficits. 
