
test_that("ReturnPeriodBiGPD works correctly", {
  expect_equal(ReturnPeriodBiGPD(c(0.9, 0.9, 0.001), c(1, 1, 1), 1, 100), 10)
  expect_equal(ReturnPeriodBiGPD(c(0.9, 0.9, 0.001), c(0.2, 0.2, 0.2), 5, 100), 50)
  expect_equal(ReturnPeriodBiGPD(c(1, 0.9, 0.001), c(0.2, 0.2, 0.2), 5, 100), Inf)
  expect_equal(ReturnPeriodBiGPD(c(0.9, 1, 0.001), c(0.2, 0.2, 0.2), 5, 100), Inf)
  expect_equal(ReturnPeriodBiGPD(c(0.9, 0.9, 0), c(0.2, 0.2, 0.2), 5, 100), Inf)
  expect_equal(ReturnPeriodBiGPD(c(0.99, 0.99, 0.0005), c(0.2, 0.2, 0.2), 10, 100), 84.7278119)
})

test_that("UnivariateExtremalIndex works correctly", {
  expect_equal(UnivariateExtremalIndex(c(1, 12, 4, 2, 25), 0.9999, 1, 3), 1)
})


data <- data.frame("Seine" = Seine_API_ERA5_1992_2021, "Loire" = Loire_API_ERA5_1992_2021)

BiGPDApproach(data, c(max(data$Seine), max(data$Loire)), c(4, 4), c(0.5, 1.5, 4.5, 6, -0.1), c(0.5, 1.5, 4.5, 6, -0.1), 0.95, 61, 30, 4)


data2 <- data.frame("TP" = Total_Precipitation_ERA5_1992_2021, "API" = API_ERA5_1992_2021)

BiGPDApproach(data2, c(max(data2$TP), 64.4), c(1, 4), c(2, 3, 0.1), c(0.5, 5, 20, 10, 0.1), 0.95, 92, 30, 7)
