# Clear the workspace
rm(list = ls())
library(dbhydroR)

station_ids <- c('S191_S',	'S65E_S','S65EX1_S',	'S84_S',	'S154_C',	'S71_S',	'S72_S',
                 'FISHP',	'S308.DS',	'L8.441',	'S133_P',	'S127_C',	'S127_P',	
                 'S129_C','S135_C',	'S2_P',	'S3_P','S4_P',	'S351_S',
                 'S352_S',	'S354_S', 'S129 PMP_P','S135 PMP_P')
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "SW", param = "FLOW",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "SW", param = "FLOW",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

#Paste dbkey and assign it a variable k
dbkeys <- c("91370", "91371", "91373", "91377", "91379", "91401", "91429", 
            "91473", "91508", "91510", "91513", "91599", "91608", "91656",
            "91668", "91675", "91687","15627", "15640", "15626","15642","15638")
for (i in dbkeys) {
  # Wrap the get_hydro() function call in a tryCatch block
  tryCatch({
    # Retrieve data for the dbkey
    data <- get_hydro(dbkey = i, date_min = "1990-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
    
    # Check if data is empty or contains only the "date" column
    if (ncol(data) <= 1) {
      cat("No data found for dbkey", i, "Skipping to the next dbkey.\n")
      next  # Skip to the next iteration of the loop
    }
    
    # Multiply all columns except "date" column by 0.0283168466 * 86400 to convert Flow rate from cfs to m³/day
    data[, -1] <- data[, -1] * (0.0283168466 * 86400)
    
    # Extract the column names excluding the date column
    column_names <- names(data)[-1]
    
    
    # Generate the filename based on the column names
    filename <- paste0( gsub(" ", "_", sub("_[^_]*$", "", paste(column_names, collapse = "_"))), "_cmd.csv")
    # Save data to a CSV file
    write.csv(data, file = filename)
    
    # Print a message indicating the file has been saved
    cat("CSV file", filename, "has been saved.\n")
    
    # Add a delay between requests
    Sys.sleep(10)  # Wait for 10 seconds before the next iteration
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    cat("No data found for dbkey", i, "Skipping to the next dbkey.\n")
  })
}

#S65E_Total
S65E_total = get_hydro(dbkey = c("91656", "AL760"), date_min = "1972-01-01", date_max = "2023-06-30") 
S65E_total[, -1] <- S65E_total[, -1] * (0.0283168466 * 86400)
write.csv(S65E_total,file ='C:/Work/Research/Data Analysis/Lake_O_Flow_data/Inflows_Dec22/S65E_total.csv')


