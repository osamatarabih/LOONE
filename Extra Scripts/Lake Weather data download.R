#Important note: run one by one,
#run all at one time gives error like Error in dbh_GET(servfull, query = qy) : Too Many Requests (RFC 6585) (HTTP 429).
# Clear the workspace
rm(list = ls())
library(dbhydroR)
library(dplyr)

station_ids <- c("L001", "L005", "L006", "LZ40")
#Rain
dbkeys <- get_dbkey(stationid = station_ids,  category = "WEATHER", param = "RAIN",stat = "SUM",recorder="CR10",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "SUM", category = "WEATHER", recorder="CR10",param = "RAIN",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

dbkeys <- c("16021", "12515", "12524", "13081")

for (i in dbkeys) {
  # Retrieve data for the dbkey
  data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
  
  # Extract the column names excluding the date column
  column_names <- names(data)[-1]
  
  # Generate the filename based on the column names
  filename <- paste0( paste(column_names, collapse = "_"), ".csv")
  
  # Save data to a CSV file
  write.csv(data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
  
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}

L001_RAIN_Inches <- read.csv("L001_RAIN_Inches.csv", colClasses = c("NULL", "character", "numeric"))
L005_RAIN_Inches = read.csv("L005_RAIN_Inches.csv", colClasses = c("NULL", "character", "numeric"))
L006_RAIN_Inches = read.csv("L006_RAIN_Inches.csv", colClasses = c("NULL", "character", "numeric"))
LZ40_RAIN_Inches = read.csv("LZ40_RAIN_Inches.csv", colClasses = c("NULL", "character", "numeric"))
#Replace NA values with zero
L001_RAIN_Inches[is.na(L001_RAIN_Inches)] <- 0
L005_RAIN_Inches[is.na(L005_RAIN_Inches)] <- 0
L006_RAIN_Inches[is.na(L006_RAIN_Inches)] <- 0
LZ40_RAIN_Inches[is.na(LZ40_RAIN_Inches)] <- 0
# Merge the files by the "date" column
merged_data <- merge(L001_RAIN_Inches, L005_RAIN_Inches, by = "date",all = TRUE)
merged_data <- merge(merged_data, L006_RAIN_Inches, by = "date",all = TRUE)
merged_data <- merge(merged_data, LZ40_RAIN_Inches, by = "date",all = TRUE)
# Calculate the average rainfall per day
merged_data$average_rainfall <- rowMeans(merged_data[, -1],na.rm = TRUE)

# View the updated merged data
head(merged_data)
# Save merged data as a CSV file
write.csv(merged_data, "LAKE_RAINFALL_DATA.csv", row.names = TRUE)



#ETPI
dbkeys <- get_dbkey(stationid = station_ids,  category = "WEATHER", param = "ETPI",stat = "SUM",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "SUM", category = "WEATHER", param = "ETPI",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

dbkeys <- c("UT736", "VM675", "UT743", "UT748")

for (i in dbkeys) {
  # Retrieve data for the dbkey
  data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
  
  # Extract the column names excluding the date column
  column_names <- names(data)[-1]
  
  # Generate the filename based on the column names
  filename <- paste0( paste(column_names, collapse = "_"), ".csv")
  
  # Save data to a CSV file
  write.csv(data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
  
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}

L001_ETPI_Inches <- read.csv("L001_ETPI_Inches.csv", colClasses = c("NULL", "character", "numeric"))
L005_ETPI_Inches = read.csv("L005_ETPI_Inches.csv", colClasses = c("NULL", "character", "numeric"))
L006_ETPI_Inches = read.csv("L006_ETPI_Inches.csv", colClasses = c("NULL", "character", "numeric"))
LZ40_ETPI_Inches = read.csv("LZ40_ETPI_Inches.csv", colClasses = c("NULL", "character", "numeric"))

# Replace NA values with zero
L001_ETPI_Inches[is.na(L001_ETPI_Inches)] <- 0
L005_ETPI_Inches[is.na(L005_ETPI_Inches)] <- 0
L006_ETPI_Inches[is.na(L006_ETPI_Inches)] <- 0
LZ40_ETPI_Inches[is.na(LZ40_ETPI_Inches)] <- 0
# Merge the files by the "date" column
merged_data <- merge(L001_ETPI_Inches, L005_ETPI_Inches, by = "date",all = TRUE)
merged_data <- merge(merged_data, L006_ETPI_Inches, by = "date",all = TRUE)
merged_data <- merge(merged_data, LZ40_ETPI_Inches, by = "date",all = TRUE)
# Calculate the average rainfall per day
merged_data$average_ETPI <- rowMeans(merged_data[, -1],na.rm = TRUE)

# View the updated merged data
head(merged_data)
# Save merged data as a CSV file
write.csv(merged_data, "LOONE_AVERAGE_ETPI_DATA.csv", row.names = TRUE)





#H2OT
dbkeys <- get_dbkey(stationid = station_ids,  category = "WQ", param = "H2OT",stat = "MEAN",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, category = "WQ", param = "H2OT",stat = "MEAN",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

dbkeys <- c("16031","12518","12527","16267")

for (i in dbkeys) {
  # Retrieve data for the dbkey
  data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
  
  # Extract the column names excluding the date column
  column_names <- names(data)[-1]
  
  # Generate the filename based on the column names
  filename <- paste0( paste(column_names, collapse = "_"), ".csv")
  
  # Save data to a CSV file
  write.csv(data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
  
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}








#RADP
dbkeys <- get_dbkey(stationid = station_ids,  category = "WEATHER", param = "RADP",stat = "MEAN",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "WEATHER", param = "RADP",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

dbkeys <- c("16025", "12516", "12525", "15649")

for (i in dbkeys) {
  # Retrieve data for the dbkey
  data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
  
  # Extract the column names excluding the date column
  column_names <- names(data)[-1]
  
  # Generate the filename based on the column names
  filename <- paste0( gsub(" ", "_", sub("_[^_]*$", "", paste(column_names, collapse = "_"))), ".csv")
  
  # Save data to a CSV file
  write.csv(data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
  
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}





#RADT
dbkeys <- get_dbkey(stationid = station_ids,  category = "WEATHER", param = "RADT",stat = "MEAN",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "WEATHER", param = "RADT",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

dbkeys <- c("16024", "12512", "12522", "13080")

for (i in dbkeys) {
  # Retrieve data for the dbkey
  data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
  
  # Extract the column names excluding the date column
  column_names <- names(data)[-1]
  
  # Generate the filename based on the column names
  filename <- paste0( gsub(" ", "_", sub("_[^_]*$", "", paste(column_names, collapse = "_"))), ".csv")
  
  # Save data to a CSV file
  write.csv(data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
  
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}







#AIRT
dbkeys <- get_dbkey(stationid = station_ids,  category = "WEATHER", param = "AIRT",stat = "MEAN",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "WEATHER", param = "AIRT",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

dbkeys <- c("16027", "12514", "12911", "13078")

for (i in dbkeys) {
  # Retrieve data for the dbkey
  data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
  
  # Extract the column names excluding the date column
  column_names <- names(data)[-1]
  
  # Generate the filename based on the column names
  filename <- paste0( paste(column_names, collapse = "_"), ".csv")
  
  # Save data to a CSV file
  write.csv(data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
  
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}





#WNDS
dbkeys <- get_dbkey(stationid = station_ids,  category = "WEATHER", param = "WNDS",stat = "MEAN",freq = "DA", detail.level = "full")
print(dbkeys)
dbkeys <- get_dbkey(stationid = station_ids, stat = "MEAN", category = "WEATHER", param = "WNDS",freq = "DA", detail.level = "dbkey")
print(dbkeys)
#using dbkeys to manual input and delete dbkey starting with zero because there is no recorder for those
cat(paste('"', dbkeys, '"', sep = "", collapse = ", "))

dbkeys <- c("16023", "12510", "12520", "13076")

for (i in dbkeys) {
  # Retrieve data for the dbkey
  data <- get_hydro(dbkey = i, date_min = "2000-01-01", date_max = format(Sys.Date(), "%Y-%m-%d"))
  
  # Extract the column names excluding the date column
  column_names <- names(data)[-1]
  
  # Generate the filename based on the column names
  filename <- paste0( paste(column_names, collapse = "_"), ".csv")
  
  # Save data to a CSV file
  write.csv(data, file = filename)
  
  # Print a message indicating the file has been saved
  cat("CSV file", filename, "has been saved.\n")
  
  # Add a delay between requests
  Sys.sleep(2) # Wait for 2 seconds before the next iteration
}



