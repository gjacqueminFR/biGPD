
data <- data.frame("Seine" = Seine_API_ERA5_1992_2021, "Loire" = Loire_API_ERA5_1992_2021)

CopulaApproach(data, c(max(data$Seine), max(data$Loire)), 0.95, 61, 30, 4)

CopulaApproach(data, c(max(data$Seine), max(data$Loire)), 0.95, 61, 30, 4, c(1, 1, 1)) # No declustering
#CopulaApproach(data, c(max(data$Seine), max(data$Loire)), 0.95, 61, 30, 4, c(1, 1, 0)) # No declustering only for the margins

UnivariateDeclustering(data$Seine, 61, 30, quantile(data$Seine, 0.95), 3)
