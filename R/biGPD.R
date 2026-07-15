
#' Univariate extremal index
#'
#' @param data A vector of data for which their extremal index will be estimated.
#' @param probaQuantile Data above the threshold of this probability are used to estimate the extremal index. The value has an impact over the estimation. It should be close to the maximal value, the 95th quantile is a classical value.
#' @param nbYears The number of distinct years. Default value is 1.
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @importFrom exdex dgaps
#' @return The univariate extremal index.
#' @export

UnivariateExtremalIndex <- function(data, probaQuantile, nbYears, Dparam) {
  if (missing(nbYears)) {
    nbYears <- 1
  }
  if (missing(Dparam)) {
    Dparam <- 3
  }

  threshold <- quantile(data, probaQuantile)

  if (max(data) < threshold) {
    return(1)
  } else {
    DgapsRes <- exdex::dgaps(matrix(data, ncol = nbYears), threshold, D = Dparam)
    extremalIndex <- DgapsRes$theta
    seextremalIndex <- DgapsRes$se

    if (extremalIndex + 2 * seextremalIndex > 1) {
      extremalIndex <- 1
    }

    return(extremalIndex)
  }
}


#' Bivariate extremal index
#'
#' @param frechetData1 A vector of data for which the extremal index will be estimated.
#' @param frechetData2 A vector of data for which the extremal index will be estimated.
#' @param probaQuantile Data above the threshold of this probability are used to estimate the extremal index. The value has an impact over the estimation. It should be close to 1, 0.95 is a classical value.
#' @param probas A vector of probability at which the bivariate extermal index is estimated.
#' @param nbYears The number of distinct years. Default value is 1.
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @return The univariate extremal index.
#' @export

BivariateExtremalIndex <- function(frechetData1, frechetData2, probaQuantile, probas, nbYears, Dparam) {
  maxFrechet <- pmax(- log(probas[1]) * frechetData1, - log(probas[2]) * frechetData2)
  extremalIndex <- UnivariateExtremalIndex(maxFrechet, probaQuantile, nbYears, Dparam)
  return(extremalIndex)
}


#' Estimate the parameters of the EGPD with maximul likelihood (mle).
#'
#' @param data A vector of data with the zeros filtered.
#' @param EGPDtype 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param initParam A vector of parameters to initialize the estimation. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @return The parameters of the EGPD in the following order: for type 1, (Kappa, Sigma, Xi) and for type 4, (Prob, Kappa, Delta, Sigma, Xi).
#' @export

EstimateEGPDParameters <- function(data, EGPDtype, initParam) {

  if (EGPDtype == 4) {
    EGPDfit <- fit.extgp(data, model = 4, method = "mle", init = initParam, confint = FALSE, R = 5, plots = FALSE)
    EGPDparam <- c(EGPDfit$fit$mle["prob"][[1]], EGPDfit$fit$mle["kappa"][[1]], EGPDfit$fit$mle["delta"][[1]], EGPDfit$fit$mle["sigma"][[1]], EGPDfit$fit$mle["xi"][[1]]) # (Prob, Kappa, Delta, Sigma, Xi)
  } else if (EGPDtype == 1) {
    invisible(capture.output(EGPDfit <- fit.extgp(data, model = 1, method = "mle", init = initParam, confint = FALSE, R = 5, plots = FALSE)))
    EGPDparam <- c(EGPDfit$fit$mle["kappa"][[1]], EGPDfit$fit$mle["sigma"][[1]], EGPDfit$fit$mle["xi"][[1]]) # (Kappa, Sigma, Xi)
  }

  return(EGPDparam)
}


#' Compute the probability to be under the return level value with an EGPD.
#'
#' @param returnLevel The value at which the cdf is computed.
#' @param EGPDtype 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param EGPDParams The EGPD parameters. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaZero The probability to be exactly zero.
#' @return The cdf value at the return level value for an EGPD.
#' @export

ProbaEGPD <- function(returnLevel, EGPDtype, EGPDParams, probaZero) {

  if (EGPDtype == 4) {
    proba <- probaZero + (1 - probaZero) * mev::pextgp(returnLevel, type = 4, prob = EGPDParams[1], kappa = EGPDParams[2], delta = EGPDParams[3], sigma = EGPDParams[4], xi = EGPDParams[5])
  } else if (EGPDtype == 1) {
    proba <- probaZero + (1 - probaZero) * mev::pextgp(returnLevel, type = 1, kappa = EGPDParams[1], sigma = EGPDParams[2], xi = EGPDParams[3])
  }

  return(proba)
}


#' Compute the value corresponding to the cdf probability with an EGPD.
#'
#' @param proba The probability at which the inverse cdf is computed.
#' @param EGPDtype 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param EGPDParams The EGPD parameters. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaZero The probability to be exactly zero.
#' @return The value at which the cdf equals proba for an EGPD.
#' @export

ValueEGPD <- function(proba, EGPDtype, EGPDParams, probaZero) {

  if (EGPDtype == 4) {
    returnLevel <- mev::qextgp((proba - probaZero) / (1 - probaZero), type = 4, prob = EGPDParams[1], kappa = EGPDParams[2], delta = EGPDParams[3], sigma = EGPDParams[4], xi = EGPDParams[5])
  } else if (EGPDtype == 1) {
    returnLevel <- mev::qextgp((proba - probaZero) / (1 - probaZero), type = 1, kappa = EGPDParams[1], sigma = EGPDParams[2], xi = EGPDParams[3])
  }

  return(returnLevel)
}


#' Transform data assumed to follow an EGPD distribution to a unit exponential distribution.
#'
#' @param data A vector of data with the zeros filtered.
#' @param EGPDtype 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param EGPDParams The EGPD parameters. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaZero The probability to be exactly zero.
#' @return A data vector.
#' @export

EGPDtoExp <- function(data, EGPDtype, EGPDParams, probaZero) {

  if (EGPDtype == 4) {
    expData <- -log(1 - (probaZero + (1 - probaZero) * mev::pextgp(data, type = 4, prob = EGPDParams[1], kappa = EGPDParams[2], delta = EGPDParams[3], sigma = EGPDParams[4], xi = EGPDParams[5])))
  } else if (EGPDtype == 1) {
    expData <- -log(1 - (probaZero + (1 - probaZero) * mev::pextgp(data, type = 1, kappa = EGPDParams[1], sigma = EGPDParams[2], xi = EGPDParams[3])))
  }

  return(expData)
}


#' Transform data assumed to follow an EGPD distribution to a unit Fréchet distribution.
#'
#' @param data A vector of data with the zeros filtered.
#' @param EGPDtype 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param EGPDParams The EGPD parameters. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaZero The probability to be exactly zero.
#' @return A data vector.
#' @export

EGPDtoFrechet <- function(data, EGPDtype, EGPDParams, probaZero) {

  if (EGPDtype == 4) {
    frechetData <- - 1 / log(probaZero + (1 - probaZero) * mev::pextgp(data, type = 4, prob = EGPDParams[1], kappa = EGPDParams[2], delta = EGPDParams[3], sigma = EGPDParams[4], xi = EGPDParams[5]))
  } else if (EGPDtype == 1) {
    frechetData <- - 1 / log(probaZero + (1 - probaZero) * mev::pextgp(data, type = 1, kappa = EGPDParams[1], sigma = EGPDParams[2], xi = EGPDParams[3]))
  }

  return(frechetData)
}


#' Bivariate return period
#'
#' @param probas A vector of size 3, with the 2 univariate probabilities first and then the probability of bivariate exceedence. The first two probabilities are likely close to 1, whereas the third one is likely close to 0.
#' @param extremalIndexes A vector of size 3, with the 2 univariate extremal indexes first and then the bivariate extremal index.
#' @param h The parameter of non-concurrence (integer).
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @return The bivariate return period.
#' @export

ReturnPeriodBiGPD <- function(probas, extremalIndexes, h, nbDaysPerYear, probaOccurrence) {
  if (missing(probaOccurrence)) {
    probaOccurrence <- 1 - exp(-1)
  }

  if (probas[1] == 1 || probas[2] == 1 || probas[3] == 0) {
    returnPeriod <- Inf
  } else {
    returnPeriod <- -log(1 - probaOccurrence) * h / (nbDaysPerYear * (1 - probas[1]^(extremalIndexes[1] * h) - probas[2]^(extremalIndexes[2] * h) + (probas[1] + probas[2] - 1 + probas[3])^(extremalIndexes[3] * h)))
  }
  return(returnPeriod)
}


#' Compute the bivariate excess probability Fbar(x1,x2).
#'
#' @param probas A vector of size 2, with the 2 univariate probabilities. The two probabilities are likely close to 1.
#' @param empCDF The empirical cdf of Delta.
#' @param probaQuantile Data above the thresholds of this probability are used to compute the probability. The value has a small impact over the probability. It should be close to 1, 0.95 is a classical value. The same value is used for both margins for simplicity.
#' @param FbarU1U2 The probability to be above both thresholds (quantiles of probability probaQuantile).
#' @importFrom cubature cubintegrate
#' @return The bivariate excess probability.
#' @export

BivariateExceedenceProbability <- function(probas, empCDF, probaQuantile, FbarU1U2) {

  IntegrandPositive <- function(t) {
    exp(-t) * empCDF(t)
  }
  IntegrandNegative <- function(t) {
    exp(t) * empCDF(t)
  }
  Integrand <- function(t) {
    exp(-t) * (empCDF(t) - empCDF(-t))
  }

  ExpValue1 <- - log(1 - probas[1])
  ExpValue2 <- - log(1 - probas[2])

  integralU1U2 <- cubature::cubintegrate(Integrand, lower = 0, upper = Inf, method = "pcubature")$integral


  errorIntegral <- FALSE
  tryCatch({
    integralX1X2 <- ((1 - probas[2]) * cubature::cubintegrate(IntegrandPositive, lower = max(0, ExpValue1 - ExpValue2), upper = Inf, method = "pcubature")$integral - (1 - probas[1]) * cubature::cubintegrate(IntegrandNegative, lower = - Inf, upper = min(0, ExpValue1 - ExpValue2), method = "pcubature")$integral)
  }, error = function(e) {errorIntegral <<- TRUE})

  if (errorIntegral) {
    Fbar <- 0
  } else {
    Fbar <- integralX1X2 * FbarU1U2 / (integralU1U2 * (1 - probaQuantile))
  }

  return(Fbar)
}


#' Compute the bivariate excess probability Fbar(x1,x2).
#'
#' @param probas A vector of size 2, with the 2 univariate probabilities. The two probabilities are likely close to 1.
#' @param empCDF The empirical cdf of Delta.
#' @param probaQuantile Data above the thresholds of this probability are used to compute the probability. The value has a small impact over the probability. It should be close to 1, 0.95 is a classical value. The same value is used for both margins for simplicity.
#' @param FU1U2 The probability to be below both thresholds (quantiles of probability probaQuantile).
#' @importFrom cubature cubintegrate
#' @return The bivariate excess probability.
#' @export

BivariateExceedenceProbability2 <- function(probas, empCDF, probaQuantile, FU1U2) {

  IntegrandPositive <- function(t) {
    exp(-t) * empCDF(t)
  }
  IntegrandNegative <- function(t) {
    exp(t) * empCDF(t)
  }
  Integrand <- function(t) {
    exp(-t) * (empCDF(t) - empCDF(-t))
  }

  ExpValueTP <- - log(1 - probas[1])
  ExpValueAPI <- - log(1 - probas[2])

  integralU1U2 <- cubature::cubintegrate(Integrand, lower = 0, upper = Inf, method = "pcubature")$integral
  print(integralU1U2)
  print((FU1U2 + 2 * (1 - probaQuantile) - 1) / (1 - FU1U2))
  print((FU1U2 + 2 * (1 - probaQuantile) - 1) / (1 - FU1U2) - integralU1U2)

  integralX1X2 <- ((1 - probas[2]) * cubature::cubintegrate(IntegrandPositive, lower = max(0, ExpValueTP - ExpValueAPI), upper = Inf, method = "pcubature")$integral - (1 - probas[1]) * cubature::cubintegrate(IntegrandNegative, lower = - Inf, upper = min(0, ExpValueTP - ExpValueAPI), method = "pcubature")$integral)

  Fbar <- integralX1X2 * (1 - FU1U2) / (1 - probaQuantile)
  print("here")

  return(Fbar)
}


#' Example of biGPD approach. Compute univariate and bivariate return periods from 2 time series.
#'
#' @param data A dataframe with two columns, one for each time series.
#' @param returnLevels A vector of size 2, one for each time series. Return periods correspond to the return levels.
#' @param EGPDtypes A vector of size 2, one for each time series. Values are 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param initParams1 A vector of parameters to initialize the estimation for the first time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param initParams2 A vector of parameters to initialize the estimation for the second time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaQuantile A float between 0 and 1. Data above the thresholds of this probability are used to compute the bivariate exceedence probability. The value has a small impact over the probability. It should be close to 1, 0.95 is a classical value. The same value is used for both margins for simplicity.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years. Default value is 1.
#' @param h The parameter of non-concurrence (integer).
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @return A list with, in that order: the return period of first variable, the return period of the second variable, the bivariate return period, the non-concurrent joint excess probability, chi and chiBarre.
#' @export

BiGPDApproach <- function(data, returnLevels, EGPDtypes, initParams1, initParams2, probaQuantile, nbDaysPerYear, nbYears, h, Dparam, probaOccurrence) {
  print("start1")
  if (missing(probaOccurrence)) {
    probaOccurrence <- 1 - exp(-1)
  }

  threshold1 <- quantile(data[, 1], probaQuantile)
  threshold2 <- quantile(data[, 2], probaQuantile)

  extremalIndex1 <- UnivariateExtremalIndex(data[, 1], probaQuantile, nbYears, Dparam)
  extremalIndex2 <- UnivariateExtremalIndex(data[, 2], probaQuantile, nbYears, Dparam)
  print("Univariate Extremal Index OK")

  # Get the probability of no rain, and filter the zeros
  probaZero1 <- length(data[, 1][data[, 1] == 0]) / length(data[, 1])
  dataFiltered1 <- data[, 1][data[, 1] > 0]
  probaZero2 <- length(data[, 2][data[, 2] == 0]) / length(data[, 2])
  dataFiltered2 <- data[, 2][data[, 2] > 0]

  EGPDparam1 <- EstimateEGPDParameters(dataFiltered1, EGPDtypes[1], initParams1)
  EGPDparam2 <- EstimateEGPDParameters(dataFiltered2, EGPDtypes[2], initParams2)
  print("EGPD parameters OK")

  # Calculate univariate return periods
  proba1 <- ProbaEGPD(returnLevels[1], EGPDtypes[1], EGPDparam1, probaZero1)
  proba2 <- ProbaEGPD(returnLevels[2], EGPDtypes[2], EGPDparam2, probaZero2)

  returnPeriod1 <- -log(1 - probaOccurrence) / (nbDaysPerYear * extremalIndex1 * (1 - proba1))
  returnPeriod2 <- -log(1 - probaOccurrence) / (nbDaysPerYear * extremalIndex2 * (1 - proba2))
  print("Univariate return period OK")

  # Constant estimated with frequency
  dataAbove <- data[data[, 1] >= threshold1 & data[, 2] >= threshold2, ]
  FbarU1U2 <- length(dataAbove[, 1]) / length(data[, 1])

  if (returnPeriod1 == Inf || returnPeriod2 == Inf) {
    returnPeriodBiv <- Inf
    HbarX1X2 <- 0
  } else {

    # Transformation to exponential margins
    dataBiv <- data[data[, 1] > threshold1 | data[, 2] > threshold2, ]

    exp1 <- EGPDtoExp(dataBiv[, 1], EGPDtypes[1], EGPDparam1, probaZero1)
    exp2 <- EGPDtoExp(dataBiv[, 2], EGPDtypes[2], EGPDparam2, probaZero2)

    thresholdExp1 <- EGPDtoExp(threshold1, EGPDtypes[1], EGPDparam1, probaZero1)
    thresholdExp2 <- EGPDtoExp(threshold2, EGPDtypes[2], EGPDparam2, probaZero2)

    # ECDF
    delta <- (exp1 - thresholdExp1) - (exp2 - thresholdExp2)
    ECDFDelta <- ecdf(delta)
    print("ECDF OK")

    # Estimate bivariate extremal index
    frechet1 <- EGPDtoFrechet(data[, 1], EGPDtypes[1], EGPDparam1, probaZero1)
    frechet2 <- EGPDtoFrechet(data[, 2], EGPDtypes[2], EGPDparam2, probaZero2)

    extremalIndexBiv <- BivariateExtremalIndex(frechet1, frechet2, probaQuantile, c(proba1, proba2), nbYears, Dparam)
    print("Bivariate Extremal Index OK")

    # Calculate bivariate return periods

    if (proba1 == 1 || proba2 == 1) {
      FbarX1X2 <- 0
    } else {
      FbarX1X2 <- BivariateExceedenceProbability(c(proba1, proba2), ECDFDelta, probaQuantile, FbarU1U2)
      #FbarX1X2 <- BivariateExceedenceProbability2(c(proba1, proba2), ECDFDelta, probaQuantile, FU1U2)
    }

    # Constraint on the bivariate extremal index
    minExtremalIndexBiv <- max(extremalIndex1 * (1 - proba1) / (2 - proba1 - proba2 - FbarX1X2), extremalIndex2 * (1 - proba2) / (2 - proba1 - proba2 - FbarX1X2))
    maxExtremalIndexBiv <- extremalIndex1 * (1 - proba1) / (2 - proba1 - proba2 - FbarX1X2) + extremalIndex2 * (1 - proba2) / (2 - proba1 - proba2 - FbarX1X2)
    if (extremalIndexBiv < minExtremalIndexBiv) {
      extremalIndexBiv <- minExtremalIndexBiv
    } else if (extremalIndexBiv > maxExtremalIndexBiv) {
      extremalIndexBiv <- maxExtremalIndexBiv
    }

    returnPeriodBiv <- ReturnPeriodBiGPD(c(proba1, proba2, FbarX1X2), c(extremalIndex1, extremalIndex2, extremalIndexBiv), h, nbDaysPerYear, probaOccurrence)
    print("Bivariate return period OK")

    # Probability of non-concurrent excess
    HbarX1X2 <- h / (nbDaysPerYear * returnPeriodBiv)
  }

  chi <- FbarU1U2 / (1 - probaQuantile)
  #chi <- (FU1U2 + 2 * (1 - probaQuantile) - 1) / (1 - probaQuantile)
  chiBarre <- (log(1 - 0.9999) - log(chi)) / (log(1 - 0.9999) + log(chi))

  return(c(returnPeriod1, returnPeriod2, returnPeriodBiv, HbarX1X2, chi, chiBarre))
}


#' Bivariate return period when searching for return levels.
#'
#' @param probaRL A float between 0 and 1. The return levels correspond to the quantiles of this probability.
#' @param probaQuantile A float between 0 and 1. Data above the thresholds of this probability are used to compute the bivariate exceedence probability. The value has a small impact over the probability. It should be close to 1, 0.95 is a classical value. The same value is used for both margins for simplicity.
#' @param FbarU1U2 The probability to be above both thresholds (quantiles of probability probaQuantile).
#' @param extremalIndex1 The univariate extremal indexe of the first variable.
#' @param extremalIndex2 The univariate extremal indexe of the second variable.
#' @param extremalIndexBivOG The bivariate extremal index.
#' @param h The parameter of non-concurrence (integer).
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @return The bivariate return period.
#' @export

ReturnPeriodBiGPDReturnLevels <- function(probaRL, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) {
  FbarX1X2 <- FbarU1U2 * (1 - probaRL) / (1 - probaQuantile)

  # Constraint on the bivariate extremal index
  minExtremalIndexBiv <- max(extremalIndex1 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2), extremalIndex2 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2))
  maxExtremalIndexBiv <- extremalIndex1 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2) + extremalIndex2 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2)
  if (extremalIndexBivOG < minExtremalIndexBiv) {
    extremalIndexBiv <- minExtremalIndexBiv
  } else if (extremalIndexBivOG > maxExtremalIndexBiv) {
    extremalIndexBiv <- maxExtremalIndexBiv
  } else {
    extremalIndexBiv <- extremalIndexBivOG
  }

  returnPeriodBiv <- ReturnPeriodBiGPD(c(probaRL, probaRL, FbarX1X2), c(extremalIndex1, extremalIndex2, extremalIndexBiv), h, nbDaysPerYear, probaOccurrence)

  return(returnPeriodBiv)
}

#' Example of biGPD approach to get return levels from return period. Compute the two return levels with equal univariate probability from a given bivariate return period..
#'
#' @param data A dataframe with two columns, one for each time series.
#' @param returnPeriod An integer, in years.
#' @param EGPDtypes A vector of size 2, one for each time series. Values are 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param initParams1 A vector of parameters to initialize the estimation for the first time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param initParams2 A vector of parameters to initialize the estimation for the second time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaQuantile A float between 0 and 1. Data above the thresholds of this probability are used to compute the bivariate exceedence probability. The value has a small impact over the probability. It should be close to 1, 0.95 is a classical value. The same value is used for both margins for simplicity.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years. Default value is 1.
#' @param h The parameter of non-concurrence (integer).
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @return A list with, in that order: the return level of first variable, the return level of the second variable, chi and chiBarre.
#' @export

BiGPDApproachReturnLevels <- function(data, returnPeriod, EGPDtypes, initParams1, initParams2, probaQuantile, nbDaysPerYear, nbYears, h, Dparam, probaOccurrence) {
  print("start1")
  if (missing(probaOccurrence)) {
    probaOccurrence <- 1 - exp(-1)
  }

  threshold1 <- quantile(data[, 1], probaQuantile)
  threshold2 <- quantile(data[, 2], probaQuantile)

  extremalIndex1 <- UnivariateExtremalIndex(data[, 1], probaQuantile, nbYears, Dparam)
  extremalIndex2 <- UnivariateExtremalIndex(data[, 2], probaQuantile, nbYears, Dparam)
  print("Univariate Extremal Index OK")

  # Get the probability of no rain, and filter the zeros
  probaZero1 <- length(data[, 1][data[, 1] == 0]) / length(data[, 1])
  dataFiltered1 <- data[, 1][data[, 1] > 0]
  probaZero2 <- length(data[, 2][data[, 2] == 0]) / length(data[, 2])
  dataFiltered2 <- data[, 2][data[, 2] > 0]

  EGPDparam1 <- EstimateEGPDParameters(dataFiltered1, EGPDtypes[1], initParams1)
  EGPDparam2 <- EstimateEGPDParameters(dataFiltered2, EGPDtypes[2], initParams2)
  print("EGPD parameters OK")

  # Transformation to exponential margins
  dataBiv <- data[data[, 1] > threshold1 | data[, 2] > threshold2, ]

  exp1 <- EGPDtoExp(dataBiv[, 1], EGPDtypes[1], EGPDparam1, probaZero1)
  exp2 <- EGPDtoExp(dataBiv[, 2], EGPDtypes[2], EGPDparam2, probaZero2)

  thresholdExp1 <- EGPDtoExp(threshold1, EGPDtypes[1], EGPDparam1, probaZero1)
  thresholdExp2 <- EGPDtoExp(threshold2, EGPDtypes[2], EGPDparam2, probaZero2)

  # ECDF
  delta <- (exp1 - thresholdExp1) - (exp2 - thresholdExp2)
  ECDFDelta <- ecdf(delta)
  print("ECDF OK")

  # Estimate bivariate extremal index
  frechet1 <- EGPDtoFrechet(data[, 1], EGPDtypes[1], EGPDparam1, probaZero1)
  frechet2 <- EGPDtoFrechet(data[, 2], EGPDtypes[2], EGPDparam2, probaZero2)

  extremalIndexBivOG <- BivariateExtremalIndex(frechet1, frechet2, probaQuantile, c(0.95, 0.95), nbYears, Dparam) # The value of the univariate probability is not necessary here, the only thing that matters is that proba1 = proba2.
  print("Bivariate Extremal Index OK")

  # Calculate return levels
  dataAbove <- data[data[, 1] >= threshold1 & data[, 2] >= threshold2, ]
  FbarU1U2 <- length(dataAbove[, 1]) / length(data[, 1])

  iteration <- 0
  probaRLmin <- probaQuantile
  probaRLmax <- 1

  if (ReturnPeriodBiGPDReturnLevels(probaRLmin, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) > returnPeriod) {
    probaRL <- probaRLmin
    print("Problem !!")
  } else {
    probaRL <- (probaRLmin + probaRLmax) / 2
    res <- ReturnPeriodBiGPDReturnLevels(probaRL, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) - returnPeriod

    while (abs(res) > 1 && iteration < 1000) {
      if (res > 0) {
        probaRLmax <- probaRL
        probaRL <- (probaRLmin + probaRLmax) / 2
      } else {
        probaRLmin <- probaRL
        probaRL <- (probaRLmin + probaRLmax) / 2
      }
      res <- ReturnPeriodBiGPDReturnLevels(probaRL, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) - returnPeriod
      iteration <- iteration + 1
    }
  }

  returnLevel1 <- ValueEGPD(probaRL, EGPDtypes[1], EGPDparam1, probaZero1)
  returnLevel2 <- ValueEGPD(probaRL, EGPDtypes[2], EGPDparam2, probaZero2)
  print("Return levels OK")

  # Chi and chi barre values
  chi <- FbarU1U2 / (1 - probaQuantile)
  chiBarre <- (log(1 - 0.9999) - log(chi)) / (log(1 - 0.9999) + log(chi))

  return(c(returnLevel1, returnLevel2, chi, chiBarre))
}
