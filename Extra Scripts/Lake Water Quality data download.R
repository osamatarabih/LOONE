# Clear the workspace
rm(list = ls())

# Load the required libraries
library(rio)
library(dbhydroR)

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
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
  # Add a delay between requests
  
}

#___________PHOSPHATE, ORTHO AS P_____________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("PHOSPHATE, ORTHO AS P")

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
  # Add a delay between requests
  
}

#____________"AMMONIA-N"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
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
  # Add a delay between requests
  
}

#____________"NITRATE+NITRITE-N"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
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
  # Add a delay between requests
  
}

#____________"TOTAL NITROGEN"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN HILR"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN HILR")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN HTYR"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN HTYR")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN LA"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN LA")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN LF"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN LF")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN LR"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN LR")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN LW"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN LW")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN LY"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN LY")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN RR"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN RR")

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
  # Add a delay between requests
  
}


#____________"MICROCYSTIN WR"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN WR")

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
  # Add a delay between requests
  
}

#____________"MICROCYSTIN YR"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("MICROCYSTIN YR")

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
  # Add a delay between requests
  
}

#____________"CHLOROPHYLL-A"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("CHLOROPHYLL-A")

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
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}

#____________"CHLOROPHYLL-A(LC)"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("CHLOROPHYLL-A(LC)")

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
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}


#____________"CHLOROPHYLL-A, CORRECTED"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("CHLOROPHYLL-A, CORRECTED")

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
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}

#____________"DISSOLVED OXYGEN"______________

# Specify the station IDs, date range, and test names
station_ids <- c("L001", "L004", "L005", "L006", "L007", "L008", "LZ40")
date_min <- "1950-01-01"
date_max <- format(Sys.Date(), "%Y-%m-%d")
test_names <- c("DISSOLVED OXYGEN")

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
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}

