# Clear the workspace
rm(list = ls())
library(dbhydroR)

station_ids <- c('S308.DS',	'S77_S',	'L8.441',	'S127_C',	'S129_C',	'S135_C',	
                 'S351_S',	'S352_S',	'S354_S',	'INDUST', 'S79','S80','S2_NNR','S3','S48_S','S49_S')
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "SW", recorder= "PREF",param = "FLOW",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "SW", recorder= "PREF",param = "FLOW",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

#Paste dbkey and assign it a variable k
dbkeys <- c("91370", "91373", "91379", "91508", "91510", "91513", "91677"
            ,"15628", "15640", "15626", "00865","JW224","00436","15018","91606","JW223")
for (i in dbkeys) {
  # Wrap the get_hydro() function call in a tryCatch block
  tryCatch({
    # Retrieve data for the dbkey
    data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
    
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



