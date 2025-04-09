
data <- data.frame("Seine" = Seine_API_ERA5_1992_2021, "Loire" = Loire_API_ERA5_1992_2021)

CopulaApproach(data, c(max(data$Seine), max(data$Loire)), 0.95, 61, 30, 4)

#CopulaApproach(data, c(max(data$Seine), max(data$Loire)), 0.95, 61, 30, 4, c(1, 1, 1)) # No declustering
#CopulaApproach(data, c(max(data$Seine), max(data$Loire)), 0.95, 61, 30, 4, c(1, 1, 0)) # No declustering only for the margins


data2 <- data.frame("TP" = Total_Precipitation_ERA5_1992_2021, "API" = API_ERA5_1992_2021)

CopulaApproach(data2, c(max(data2$TP), 64.38838), 0.95, 92, 30, 7, c(2, 7, 2))
