# Clear the workspace
rm(list = ls())

# Load the required libraries
library(rio)
library(dbhydroR)
#Stage Data
LO_Stage = get_hydro(dbkey = c("16022", "12509","12519","16265","15611"), date_min = "1950-01-01",
                     date_max = format(Sys.Date(), "%Y-%m-%d"))
write.csv(LO_Stage,file ='C:/Work/LOONE_Data/LO_Stage.csv')

#WCAs
Stg_3ANW = get_hydro(dbkey = "LA369", date_min = "1972-01-01", date_max = "2023-04-30")
write.csv(Stg_3ANW,file ='C:/Work/Research/Data Analysis/Lake_O_Flow_data/Inflows_Dec22/Stg_3ANW.csv')
Stg_2A17 = get_hydro(dbkey = "16531", date_min = "1972-01-01", date_max = "2023-04-30")
write.csv(Stg_2A17,file ='C:/Work/Research/Data Analysis/Lake_O_Flow_data/Inflows_Dec22/Stg_2A17.csv')
Stg_3A3 = get_hydro(dbkey = "16532", date_min = "1972-01-01", date_max = "2023-04-30")
write.csv(Stg_3A3,file ='C:/Work/Research/Data Analysis/Lake_O_Flow_data/Inflows_Dec22/Stg_3A3.csv')
Stg_3A4 = get_hydro(dbkey = "16537", date_min = "1972-01-01", date_max = "2023-04-30")
write.csv(Stg_3A4,file ='C:/Work/Research/Data Analysis/Lake_O_Flow_data/Inflows_Dec22/Stg_3A4.csv')
Stg_3A28 = get_hydro(dbkey = "16538", date_min = "1972-01-01", date_max = "2023-04-30")
write.csv(Stg_3A28,file ='C:/Work/Research/Data Analysis/Lake_O_Flow_data/Inflows_Dec22/Stg_3A28.csv')