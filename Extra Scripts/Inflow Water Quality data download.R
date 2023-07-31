# Clear the workspace
rm(list = ls())

# Load the required libraries
library(rio)
library(dbhydroR)

# Specify the station IDs, date range, and test names
station_ids <- c('S191',	'S65E',	'S84',	'S154',	'S71',	'S72',
                 'S4',	'FECSR78',	'S308C',	'CULV10A',	'S133',
                 'S127',	'S135')
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("PHOSPHATE, TOTAL AS P")

# Loop over the station IDs
for (station_id in station_ids) {
  # Retrieve water quality data for the current station ID
  water_quality_data <- get_wq(
    station_id = station_id,
    date_min = date_min,
    date_max = date_max,
    test_name = test_names
  )
  
  # Convert negative values to NA
  water_quality_data[water_quality_data < 0] <- NA
  
  # Calculate the number of days from the minimum date plus 8
  water_quality_data$days <- as.integer(difftime(water_quality_data$date, min(water_quality_data$date), units = "days")) + as.integer(format(min(water_quality_data$date), "%d"))
  
  # Generate the filename based on the station ID
  filename <- paste0("water_quality_",station_id,"_",test_names, ".csv")
  
  # Save data to a CSV file
  write.csv(water_quality_data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
}

#____________"AMMONIA-N"______________

test_names <- c("AMMONIA-N")

# Loop over the station IDs
for (station_id in station_ids) {
  # Retrieve water quality data for the current station ID
  water_quality_data <- get_wq(
    station_id = station_id,
    date_min = date_min,
    date_max = date_max,
    test_name = test_names
  )
  
  # Convert negative values to NA
  water_quality_data[water_quality_data < 0] <- NA
  
  # Calculate the number of days from the minimum date plus 8
  water_quality_data$days <- as.integer(difftime(water_quality_data$date, min(water_quality_data$date), units = "days")) + as.integer(format(min(water_quality_data$date), "%d"))
  
  # Generate the filename based on the station ID
  filename <- paste0("water_quality_",station_id,"_",test_names, ".csv")
  
  # Save data to a CSV file
  write.csv(water_quality_data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
}

#____________"NITRATE+NITRITE-N"______________

test_names <- c("NITRATE+NITRITE-N")

# Loop over the station IDs
for (station_id in station_ids) {
  # Retrieve water quality data for the current station ID
  water_quality_data <- get_wq(
    station_id = station_id,
    date_min = date_min,
    date_max = date_max,
    test_name = test_names
  )
  
  # Convert negative values to NA
  water_quality_data[water_quality_data < 0] <- NA
  
  # Calculate the number of days from the minimum date plus 8
  water_quality_data$days <- as.integer(difftime(water_quality_data$date, min(water_quality_data$date), units = "days")) + as.integer(format(min(water_quality_data$date), "%d"))
  
  # Generate the filename based on the station ID
  filename <- paste0("water_quality_",station_id,"_",test_names, ".csv")
  
  # Save data to a CSV file
  write.csv(water_quality_data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
}

#____________"TOTAL NITROGEN"______________

test_names <- c("TOTAL NITROGEN")

# Loop over the station IDs
for (station_id in station_ids) {
  # Retrieve water quality data for the current station ID
  water_quality_data <- get_wq(
    station_id = station_id,
    date_min = date_min,
    date_max = date_max,
    test_name = test_names
  )
  
  # Convert negative values to NA
  water_quality_data[water_quality_data < 0] <- NA
  
  # Calculate the number of days from the minimum date plus 8
  water_quality_data$days <- as.integer(difftime(water_quality_data$date, min(water_quality_data$date), units = "days")) + as.integer(format(min(water_quality_data$date), "%d"))
  
  # Generate the filename based on the station ID
  filename <- paste0("water_quality_",station_id,"_",test_names, ".csv")
  
  # Save data to a CSV file
  write.csv(water_quality_data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
}


#____________"CHLOROPHYLL-A"______________

# Specify the station IDs, date range, and test names
station_ids <- c('S65E',	'S84',	'S154',	'S71',	'S72',
                 'S4',	'FECSR78',	'S308C',	'CULV10A',	'S133',
                 'S127',	'S135', 'S191')
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- "CHLOROPHYLL-A"

# Loop over the station IDs
for (station_id in station_ids) {
  # Retrieve water quality data for the current station ID
  water_quality_data <- tryCatch(
    get_wq(
      station_id = station_id,
      date_min = date_min,
      date_max = date_max,
      test_name = test_names
    ),
    error = function(e) NULL
  )
  
  # Check if data is available for the current station ID and test name
  if (!is.null(water_quality_data) && nrow(water_quality_data) > 0) {
    # Convert negative values to NA
    water_quality_data[water_quality_data < 0] <- NA
    
    # Convert the vector to a data frame
    water_quality_data <- as.data.frame(water_quality_data)
    
    # Calculate the number of days from the minimum date plus 8
    water_quality_data$days <- as.integer(difftime(water_quality_data$date, min(water_quality_data$date), units = "days")) + as.integer(format(min(water_quality_data$date), "%d"))
    
    # Generate the filename based on the station ID
    filename <- paste0("water_quality_", station_id, "_", test_names, ".csv")
    
    # Save data to a CSV file
    write.csv(water_quality_data, file = filename)
    
    # Print a message indicating the file has been saved
    cat("CSV file", filename, "has been saved.\n")
  } else {
    # Print a message indicating no data was found for the current station ID and test name
    cat("No data found for station ID", station_id, "and test name", test_names, "\n")
  }
  
  Sys.sleep(1) # Wait for 1 seconds before the next iteration
}


#____________"CHLOROPHYLL-A(LC)"______________

# Specify the station IDs, date range, and test names
station_ids <- c('S65E',	'S84',	'S154',	'S71',	'S72',
                 'S4',	'FECSR78',	'S308C',	'CULV10A',	'S133',
                 'S127',	'S135', 'S191')
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- "CHLOROPHYLL-A(LC)"

# Loop over the station IDs
for (station_id in station_ids) {
  # Retrieve water quality data for the current station ID
  water_quality_data <- tryCatch(
    get_wq(
      station_id = station_id,
      date_min = date_min,
      date_max = date_max,
      test_name = test_names
    ),
    error = function(e) NULL
  )
  
  # Check if data is available for the current station ID and test name
  if (!is.null(water_quality_data) && nrow(water_quality_data) > 0) {
    # Convert negative values to NA
    water_quality_data[water_quality_data < 0] <- NA
    
    # Convert the vector to a data frame
    water_quality_data <- as.data.frame(water_quality_data)
    
    # Calculate the number of days from the minimum date plus 8
    water_quality_data$days <- as.integer(difftime(water_quality_data$date, min(water_quality_data$date), units = "days")) + as.integer(format(min(water_quality_data$date), "%d"))
    
    # Generate the filename based on the station ID
    filename <- paste0("water_quality_", station_id, "_", test_names, ".csv")
    
    # Save data to a CSV file
    write.csv(water_quality_data, file = filename)
    
    # Print a message indicating the file has been saved
    cat("CSV file", filename, "has been saved.\n")
  } else {
    # Print a message indicating no data was found for the current station ID and test name
    cat("No data found for station ID", station_id, "and test name", test_names, "\n")
  }
  
  Sys.sleep(1) # Wait for 1 seconds before the next iteration
}

#____________"CHLOROPHYLL-A, CORRECTED"______________

# Specify the station IDs, date range, and test names
station_ids <- c('S65E',	'S84',	'S154',	'S71',	'S72',
                 'S4',	'FECSR78',	'S308C',	'CULV10A',	'S133',
                 'S127',	'S135', 'S191')
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- "CHLOROPHYLL-A, CORRECTED"

# Loop over the station IDs
for (station_id in station_ids) {
  # Retrieve water quality data for the current station ID
  water_quality_data <- tryCatch(
    get_wq(
      station_id = station_id,
      date_min = date_min,
      date_max = date_max,
      test_name = test_names
    ),
    error = function(e) NULL
  )
  
  # Check if data is available for the current station ID and test name
  if (!is.null(water_quality_data) && nrow(water_quality_data) > 0) {
    # Convert negative values to NA
    water_quality_data[water_quality_data < 0] <- NA
    
    # Convert the vector to a data frame
    water_quality_data <- as.data.frame(water_quality_data)
    
    # Calculate the number of days from the minimum date plus 8
    water_quality_data$days <- as.integer(difftime(water_quality_data$date, min(water_quality_data$date), units = "days")) + as.integer(format(min(water_quality_data$date), "%d"))
    
    # Generate the filename based on the station ID
    filename <- paste0("water_quality_", station_id, "_", test_names, ".csv")
    
    # Save data to a CSV file
    write.csv(water_quality_data, file = filename)
    
    # Print a message indicating the file has been saved
    cat("CSV file", filename, "has been saved.\n")
  } else {
    # Print a message indicating no data was found for the current station ID and test name
    cat("No data found for station ID", station_id, "and test name", test_names, "\n")
  }
  
  Sys.sleep(1) # Wait for 1 seconds before the next iteration
}

