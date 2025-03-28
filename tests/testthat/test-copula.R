
data <- data.frame("Seine" = Seine_API_ERA5_1992_2021, "Loire" = Loire_API_ERA5_1992_2021)

CopulaApproach(data, 0.95, c(max(data$Seine), max(data$Loire)), 61, 30, 4)

UnivariateDeclustering(data$Seine, 61, 30, quantile(data$Seine, 0.95), 3)
